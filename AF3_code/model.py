# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
# this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

"""AlphaFold3 model."""

from collections.abc import Iterable, Mapping
import concurrent
import dataclasses
import functools
from typing import Any, TypeAlias

from absl import logging
from alphafold3 import structure
from alphafold3.common import base_config
from alphafold3.model import confidences
from alphafold3.model import feat_batch
from alphafold3.model import features
from alphafold3.model import model_config
from alphafold3.model.atom_layout import atom_layout
from alphafold3.model.components import mapping
from alphafold3.model.components import utils
from alphafold3.model.network import atom_cross_attention
from alphafold3.model.network import confidence_head
from alphafold3.model.network import diffusion_head
from alphafold3.model.network import distogram_head
from alphafold3.model.network import evoformer as evoformer_network
from alphafold3.model.network import featurization
import haiku as hk
import jax
import jax.numpy as jnp
import numpy as np


ModelResult: TypeAlias = Mapping[str, Any]
_ScalarNumberOrArray: TypeAlias = Mapping[str, float | int | np.ndarray]


@dataclasses.dataclass(frozen=True)
class InferenceResult:
  """Postprocessed model result.

  Attributes:
    predicted_structure: Predicted protein structure.
    numerical_data: Useful numerical data (scalars or arrays) to be saved at
      inference time.
    metadata: Smaller numerical data (usually scalar) to be saved as inference
      metadata.
    debug_outputs: Additional dict for debugging, e.g. raw outputs of a model
      forward pass.
    model_id: Model identifier.
  """

  predicted_structure: structure.Structure = dataclasses.field()
  numerical_data: _ScalarNumberOrArray = dataclasses.field(default_factory=dict)
  metadata: _ScalarNumberOrArray = dataclasses.field(default_factory=dict)
  debug_outputs: Mapping[str, Any] = dataclasses.field(default_factory=dict)
  model_id: bytes = b''


def get_predicted_structure(
    result: ModelResult, batch: feat_batch.Batch
) -> structure.Structure:
  """Creates the predicted structure and ion preditions.

  Args:
    result: model output in a model specific layout
    batch: model input batch

  Returns:
    Predicted structure.
  """
  model_output_coords = result['diffusion_samples']['atom_positions']

  # Rearrange model output coordinates to the flat output layout.
  model_output_to_flat = atom_layout.compute_gather_idxs(
      source_layout=batch.convert_model_output.token_atoms_layout,
      target_layout=batch.convert_model_output.flat_output_layout,
  )
  pred_flat_atom_coords = atom_layout.convert(
      gather_info=model_output_to_flat,
      arr=model_output_coords,
      layout_axes=(-3, -2),
  )

  predicted_lddt = result.get('predicted_lddt')

  if predicted_lddt is not None:
    pred_flat_b_factors = atom_layout.convert(
        gather_info=model_output_to_flat,
        arr=predicted_lddt,
        layout_axes=(-2, -1),
    )
  else:
    # Handle models which don't have predicted_lddt outputs.
    pred_flat_b_factors = np.zeros(pred_flat_atom_coords.shape[:-1])

  (missing_atoms_indices,) = np.nonzero(model_output_to_flat.gather_mask == 0)
  if missing_atoms_indices.shape[0] > 0:
    missing_atoms_flat_layout = batch.convert_model_output.flat_output_layout[
        missing_atoms_indices
    ]
    missing_atoms_uids = list(
        zip(
            missing_atoms_flat_layout.chain_id,
            missing_atoms_flat_layout.res_id,
            missing_atoms_flat_layout.res_name,
            missing_atoms_flat_layout.atom_name,
        )
    )
    logging.warning(
        'Target %s: warning: %s atoms were not predicted by the '
        'model, setting their coordinates to (0, 0, 0). '
        'Missing atoms: %s',
        batch.convert_model_output.empty_output_struc.name,
        missing_atoms_indices.shape[0],
        missing_atoms_uids,
    )

  # Put them into a structure
  pred_struc = batch.convert_model_output.empty_output_struc
  pred_struc = pred_struc.copy_and_update_atoms(
      atom_x=pred_flat_atom_coords[..., 0],
      atom_y=pred_flat_atom_coords[..., 1],
      atom_z=pred_flat_atom_coords[..., 2],
      atom_b_factor=pred_flat_b_factors,
      atom_occupancy=np.ones(pred_flat_atom_coords.shape[:-1]),  # Always 1.0.
  )
  # Set manually/differently when adding metadata.
  pred_struc = pred_struc.copy_and_update_globals(release_date=None)
  return pred_struc


def create_target_feat_embedding(
    batch: feat_batch.Batch,
    config: evoformer_network.Evoformer.Config,
    global_config: model_config.GlobalConfig,
) -> jnp.ndarray:
  """Create target feature embedding."""

  dtype = jnp.bfloat16 if global_config.bfloat16 == 'all' else jnp.float32

  with utils.bfloat16_context():
    target_feat = featurization.create_target_feat(
        batch,
        append_per_atom_features=False,
    ).astype(dtype)

    enc = atom_cross_attention.atom_cross_att_encoder(
        token_atoms_act=None,
        trunk_single_cond=None,
        trunk_pair_cond=None,
        config=config.per_atom_conditioning,
        global_config=global_config,
        batch=batch,
        name='evoformer_conditioning',
    )
    target_feat = jnp.concatenate([target_feat, enc.token_act], axis=-1).astype(
        dtype
    )

  return target_feat


def _compute_ptm(
    result: ModelResult,
    num_tokens: int,
    asym_id: np.ndarray,
    pae_single_mask: np.ndarray,
    interface: bool,
) -> np.ndarray:
  """Computes the pTM metrics from PAE."""
  return np.stack(
      [
          confidences.predicted_tm_score(
              tm_adjusted_pae=tm_adjusted_pae[:num_tokens, :num_tokens],
              asym_id=asym_id,
              pair_mask=pae_single_mask[:num_tokens, :num_tokens],
              interface=interface,
          )
          for tm_adjusted_pae in result['tmscore_adjusted_pae_global']
      ],
      axis=0,
  )


def _compute_chain_pair_iptm(
    num_tokens: int,
    asym_ids: np.ndarray,
    mask: np.ndarray,
    tm_adjusted_pae: np.ndarray,
) -> np.ndarray:
  """Computes the chain pair ipTM metrics from PAE."""
  return np.stack(
      [
          confidences.chain_pairwise_predicted_tm_scores(
              tm_adjusted_pae=sample_tm_adjusted_pae[:num_tokens],
              asym_id=asym_ids[:num_tokens],
              pair_mask=mask[:num_tokens, :num_tokens],
          )
          for sample_tm_adjusted_pae in tm_adjusted_pae
      ],
      axis=0,
  )


class Model(hk.Module):
  """Full model. Takes in data batch and returns model outputs."""

  class HeadsConfig(base_config.BaseConfig):
    diffusion: diffusion_head.DiffusionHead.Config = base_config.autocreate()
    confidence: confidence_head.ConfidenceHead.Config = base_config.autocreate()
    distogram: distogram_head.DistogramHead.Config = base_config.autocreate()

  class Config(base_config.BaseConfig):
    evoformer: evoformer_network.Evoformer.Config = base_config.autocreate()
    global_config: model_config.GlobalConfig = base_config.autocreate()
    heads: 'Model.HeadsConfig' = base_config.autocreate()
    num_recycles: int = 10
    return_embeddings: bool = False

  def __init__(self, config: Config, name: str = 'diffuser'):
    super().__init__(name=name)
    self.config = config
    self.global_config = config.global_config
    self.diffusion_module = diffusion_head.DiffusionHead(
        self.config.heads.diffusion, self.global_config
    )

  @hk.transparent
  def _sample_diffusion(
      self,
      batch: feat_batch.Batch,
      embeddings: dict[str, jnp.ndarray],
      *,
      sample_config: diffusion_head.SampleConfig,
  ) -> dict[str, jnp.ndarray]:
    denoising_step = functools.partial(
        self.diffusion_module,
        batch=batch,
        embeddings=embeddings,
        use_conditioning=True,
    )

    sample = diffusion_head.sample(
        denoising_step=denoising_step,
        batch=batch,
        key=hk.next_rng_key(),
        config=sample_config,
    )
    return sample

  def __call__(
      self, batch: features.BatchDict, key: jax.Array | None = None
  ) -> ModelResult:
    if key is None:
      key = hk.next_rng_key()

    batch = feat_batch.Batch.from_data_dict(batch)

    embedding_module = evoformer_network.Evoformer(
        self.config.evoformer, self.global_config
    )
    target_feat = create_target_feat_embedding(
        batch=batch,
        config=embedding_module.config,
        global_config=self.global_config,
    )

    def recycle_body(_, args):
      prev, key = args
      key, subkey = jax.random.split(key)
      embeddings = embedding_module(
          batch=batch,
          prev=prev,
          target_feat=target_feat,
          key=subkey,
      )
      embeddings['pair'] = embeddings['pair'].astype(jnp.float32)
      embeddings['single'] = embeddings['single'].astype(jnp.float32)
      return embeddings, key

    num_res = batch.num_res

    embeddings = {
        'pair': jnp.zeros(
            [num_res, num_res, self.config.evoformer.pair_channel],
            dtype=jnp.float32,
        ),
        'single': jnp.zeros(
            [num_res, self.config.evoformer.seq_channel], dtype=jnp.float32
        ),
        'target_feat': target_feat,
    }
    if hk.running_init():
      embeddings, _ = recycle_body(None, (embeddings, key))
    else:
      # Number of recycles is number of additional forward trunk passes.
      num_iter = self.config.num_recycles + 1
      embeddings, _ = hk.fori_loop(0, num_iter, recycle_body, (embeddings, key))

    samples = self._sample_diffusion(
        batch,
        embeddings,
        sample_config=self.config.heads.diffusion.eval,
    )

    # Compute dist_error_fn over all samples for distance error logging.
    confidence_output = mapping.sharded_map(
        lambda dense_atom_positions: confidence_head.ConfidenceHead(
            self.config.heads.confidence, self.global_config
        )(
            dense_atom_positions=dense_atom_positions,
            embeddings=embeddings,
            seq_mask=batch.token_features.mask,
            token_atoms_to_pseudo_beta=batch.pseudo_beta_info.token_atoms_to_pseudo_beta,
            asym_id=batch.token_features.asym_id,
        ),
        in_axes=0,
    )(samples['atom_positions'])

    distogram = distogram_head.DistogramHead(
        self.config.heads.distogram, self.global_config
    )(batch, embeddings)

    output = {
        'diffusion_samples': samples,
        'distogram': distogram,
        **confidence_output,
    }
    if self.config.return_embeddings:
      output['single_embeddings'] = embeddings['single']
      output['pair_embeddings'] = embeddings['pair']
      output['target_feat'] = embeddings['target_feat']
    return output

  @classmethod
  def get_inference_result(
      cls,
      batch: features.BatchDict,
      result: ModelResult,
      target_name: str = '',
  ) -> Iterable[InferenceResult]:
    """Get the predicted structure, scalars, and arrays for inference.

    This function also computes any inference-time quantities, which are not a
    part of the forward-pass, e.g. additional confidence scores. Note that this
    function is not serialized, so it should be slim if possible.

    Args:
      batch: data batch used for model inference, incl. TPU invalid types.
      result: output dict from the model's forward pass.
      target_name: target name to be saved within structure.

    Yields:
      inference_result: dataclass object that contains a predicted structure,
      important inference-time scalars and arrays, as well as a slightly trimmed
      dictionary of raw model result from the forward pass (for debugging).
    """
    del target_name
    batch = feat_batch.Batch.from_data_dict(batch)

    # Retrieve structure and construct a predicted structure.
    pred_structure = get_predicted_structure(result=result, batch=batch)

    num_tokens = batch.token_features.seq_length.item()

    pae_single_mask = np.tile(
        batch.frames.mask[:, None],
        [1, batch.frames.mask.shape[0]],
    )
    ptm = _compute_ptm(
        result=result,
        num_tokens=num_tokens,
        asym_id=batch.token_features.asym_id[:num_tokens],
        pae_single_mask=pae_single_mask,
        interface=False,
    )
    iptm = _compute_ptm(
        result=result,
        num_tokens=num_tokens,
        asym_id=batch.token_features.asym_id[:num_tokens],
        pae_single_mask=pae_single_mask,
        interface=True,
    )
    ptm_iptm_average = 0.8 * iptm + 0.2 * ptm

    asym_ids = batch.token_features.asym_id[:num_tokens]
    # Map asym IDs back to chain IDs. Asym IDs are constructed from chain IDs by
    # iterating over the chain IDs, and for each unique chain ID incrementing
    # the asym ID by 1 and mapping it to the particular chain ID. Asym IDs are
    # 1-indexed, so subtract 1 to get back to the chain ID.
    chain_ids = [pred_structure.chains[asym_id - 1] for asym_id in asym_ids]
    res_ids = batch.token_features.residue_index[:num_tokens]

    if len(np.unique(asym_ids[:num_tokens])) > 1:
      # There is more than one chain, hence interface pTM (i.e. ipTM) defined,
      # so use it.
      ranking_confidence = ptm_iptm_average
    else:
      # There is only one chain, hence ipTM=NaN, so use just pTM.
      ranking_confidence = ptm

    contact_probs = result['distogram']['contact_probs']
    # Compute PAE related summaries.
    _, chain_pair_pae_min, _ = confidences.chain_pair_pae(
        num_tokens=num_tokens,
        asym_ids=batch.token_features.asym_id,
        full_pae=result['full_pae'],
        mask=pae_single_mask,
    )
    chain_pair_pde_mean, chain_pair_pde_min = confidences.chain_pair_pde(
        num_tokens=num_tokens,
        asym_ids=batch.token_features.asym_id,
        full_pde=result['full_pde'],
    )
    intra_chain_single_pde, cross_chain_single_pde, _ = confidences.pde_single(
        num_tokens,
        batch.token_features.asym_id,
        result['full_pde'],
        contact_probs,
    )
    pae_metrics = confidences.pae_metrics(
        num_tokens=num_tokens,
        asym_ids=batch.token_features.asym_id,
        full_pae=result['full_pae'],
        mask=pae_single_mask,
        contact_probs=contact_probs,
        tm_adjusted_pae=result['tmscore_adjusted_pae_interface'],
    )
    ranking_confidence_pae = confidences.rank_metric(
        result['full_pae'],
        contact_probs * batch.frames.mask[:, None].astype(float),
    )
    chain_pair_iptm = _compute_chain_pair_iptm(
        num_tokens=num_tokens,
        asym_ids=batch.token_features.asym_id,
        mask=pae_single_mask,
        tm_adjusted_pae=result['tmscore_adjusted_pae_interface'],
    )
    # iptm_ichain is a vector of per-chain ptm values. iptm_ichain[0],
    # for example, is just the zeroth diagonal entry of the chain pair iptm
    # matrix:
    # [[x, , ],
    #  [ , , ],
    #  [ , , ]]]
    iptm_ichain = chain_pair_iptm.diagonal(axis1=-2, axis2=-1)
    # iptm_xchain is a vector of cross-chain interactions for each chain.
    # iptm_xchain[0], for example, is an average of chain 0's interactions with
    # other chains:
    # [[ ,x,x],
    #  [x, , ],
    #  [x, , ]]]
    iptm_xchain = confidences.get_iptm_xchain(chain_pair_iptm)

    predicted_distance_errors = result['average_pde']

    # Computing solvent accessible area with dssp can be slow for large
    # structures with lots of chains, so we parallelize the call.
    pred_structures = pred_structure.unstack()
    num_workers = len(pred_structures)
    with concurrent.futures.ThreadPoolExecutor(
        max_workers=num_workers
    ) as executor:
      has_clash = list(executor.map(confidences.has_clash, pred_structures))
      fraction_disordered = list(
          executor.map(confidences.fraction_disordered, pred_structures)
      )

    for idx, pred_structure in enumerate(pred_structures):
      ranking_score = confidences.get_ranking_score(
          ptm=ptm[idx],
          iptm=iptm[idx],
          fraction_disordered_=fraction_disordered[idx],
          has_clash_=has_clash[idx],
      )
      yield InferenceResult(
          predicted_structure=pred_structure,
          numerical_data={
              'full_pde': result['full_pde'][idx, :num_tokens, :num_tokens],
              'full_pae': result['full_pae'][idx, :num_tokens, :num_tokens],
              'contact_probs': contact_probs[:num_tokens, :num_tokens],
          },
          metadata={
              'predicted_distance_error': predicted_distance_errors[idx],
              'ranking_score': ranking_score,
              'fraction_disordered': fraction_disordered[idx],
              'has_clash': has_clash[idx],
              'predicted_tm_score': ptm[idx],
              'interface_predicted_tm_score': iptm[idx],
              'chain_pair_pde_mean': chain_pair_pde_mean[idx],
              'chain_pair_pde_min': chain_pair_pde_min[idx],
              'chain_pair_pae_min': chain_pair_pae_min[idx],
              'ptm': ptm[idx],
              'iptm': iptm[idx],
              'ptm_iptm_average': ptm_iptm_average[idx],
              'intra_chain_single_pde': intra_chain_single_pde[idx],
              'cross_chain_single_pde': cross_chain_single_pde[idx],
              'pae_ichain': pae_metrics['pae_ichain'][idx],
              'pae_xchain': pae_metrics['pae_xchain'][idx],
              'ranking_confidence': ranking_confidence[idx],
              'ranking_confidence_pae': ranking_confidence_pae[idx],
              'chain_pair_iptm': chain_pair_iptm[idx],
              'iptm_ichain': iptm_ichain[idx],
              'iptm_xchain': iptm_xchain[idx],
              'token_chain_ids': chain_ids,
              'token_res_ids': res_ids,
          },
          model_id=result['__identifier__'],
          debug_outputs={},
      )


class Cond_Model(Model):
  """Extend AlphaFold3 model, directly use pre-calculated embeddings without re-calculation."""

  def __init__(self, config: Model.Config, name: str = 'diffuser'):
    # Use the same name as Model to ensure parameter sharing
    super().__init__(config=config, name=name)

  def __call__(
      self, batch: features.BatchDict, key: jax.Array | None = None
  ) -> ModelResult:
    """Use pre-calculated embeddings, pocket diffusion logic, or fallback to standard AF3Model."""
    if key is None:
      key = hk.next_rng_key()

    batch_obj = feat_batch.Batch.from_data_dict(batch)
    batch_dict = batch # Original batch dict from input features.BatchDict

    # Check if pocket diffusion logic is enabled
    use_pocket_diffusion_flag = hasattr(self.config, 'use_pocket_diffusion_logic') and self.config.use_pocket_diffusion_logic
    current_embeddings = {}

    with utils.bfloat16_context():
        target_feat_base = featurization.create_target_feat(
            batch_obj, append_per_atom_features=False,
        )
        dtype = jnp.bfloat16 if self.global_config.bfloat16 == 'all' else jnp.float32
        # Note: evoformer_conditioning must be the same as in the original Model to share parameters
        enc = atom_cross_attention.atom_cross_att_encoder(
            token_atoms_act=None, trunk_single_cond=None, trunk_pair_cond=None,
            config=self.config.evoformer.per_atom_conditioning,
            global_config=self.global_config, batch=batch_obj, name='evoformer_conditioning',
        )
        calculated_target_feat = jnp.concatenate([target_feat_base, enc.token_act], axis=-1).astype(dtype)

    if use_pocket_diffusion_flag: # SiteAF3 conditional pocket diffusion logic
        print("Cond_Model: Engaging pocket sample diffusion logic.")
        if 'single_embedding' in batch_dict and 'pair_embedding' in batch_dict:
            current_embeddings['single'] = batch_dict['single_embedding'].astype(jnp.float32)
            current_embeddings['pair'] = batch_dict['pair_embedding'].astype(jnp.float32)
        else:
            print("Warning: Pocket diffusion active, but 'single_embedding' or 'pair_embedding' not in batch_dict.")

        current_embeddings['target_feat'] = calculated_target_feat
        
        if 'structure_atom_coords' not in batch_dict or \
           'structure_atom_mask' not in batch_dict or \
           'center_pos' not in batch_dict:
            raise ValueError("Cond_Model: Pocket diffusion requires 'structure_atom_coords', "
                             "'structure_atom_mask', and 'center_pos' in the input batch_dict.")
        
        ref_coords = batch_dict['structure_atom_coords']
        ref_mask = batch_dict['structure_atom_mask']
        center_position = batch_dict['center_pos']

        samples = self._pocket_sample_diffusion(
            batch_obj,
            current_embeddings, 
            reference_atom_coords=ref_coords,
            reference_atom_mask=ref_mask,
            center_pos=center_position,
            batch_dict=batch_dict,
            sample_config=self.config.heads.diffusion.eval,
        )

    elif 'single_embedding' in batch_dict and 'pair_embedding' in batch_dict: # Compatible with alphafold3
        print("Cond_Model: Using pre-calculated single/pair embeddings (standard direct use).")
        current_embeddings['single'] = batch_dict['single_embedding'].astype(jnp.float32)
        current_embeddings['pair'] = batch_dict['pair_embedding'].astype(jnp.float32)
        current_embeddings['target_feat'] = calculated_target_feat
        
        samples = self._sample_diffusion( # This is Model._sample_diffusion
            batch_obj,
            current_embeddings,
            sample_config=self.config.heads.diffusion.eval,
        )
    else:
        print("Cond_Model: No recognized conditional logic triggered, falling back to standard Model.__call__.")
        return super().__call__(batch, key) # Calls Model.__call__ which does full Evoformer etc.

    # Post-processing (confidence, distogram, output) - common to both conditional paths
    confidence_output = mapping.sharded_map(
        lambda dense_atom_positions: confidence_head.ConfidenceHead(
            self.config.heads.confidence, self.global_config
        )(
            dense_atom_positions=dense_atom_positions,
            embeddings=current_embeddings, 
            seq_mask=batch_obj.token_features.mask,
            token_atoms_to_pseudo_beta=batch_obj.pseudo_beta_info.token_atoms_to_pseudo_beta,
            asym_id=batch_obj.token_features.asym_id,
        ),
        in_axes=0,
    )(samples['atom_positions'])

    distogram = distogram_head.DistogramHead(
        self.config.heads.distogram, self.global_config
    )(batch_obj, current_embeddings)

    output = {
        'diffusion_samples': samples, 'distogram': distogram, **confidence_output,
    }
    if self.config.return_embeddings:
        if 'single' in current_embeddings: output['single_embeddings'] = current_embeddings['single']
        if 'pair' in current_embeddings: output['pair_embeddings'] = current_embeddings['pair']
        if 'target_feat' in current_embeddings: output['target_feat'] = current_embeddings['target_feat']
    return output

  @hk.transparent
  def _pocket_sample_diffusion(
      self,
      batch: feat_batch.Batch, # Standard Batch object for other AF3 components
      embeddings: dict[str, jnp.ndarray],
      reference_atom_coords: jnp.ndarray, # Explicitly passed reference coordinates
      reference_atom_mask: jnp.ndarray,   # Explicitly passed reference mask
      center_pos: jnp.ndarray,            # Added center_pos parameter
      batch_dict: dict,                  # Added batch_dict parameter
      *,
      sample_config: diffusion_head.SampleConfig,
  ) -> dict[str, jnp.ndarray]:
    '''
    Using custom diffusion, fix the receptor coordinates (from reference_atom_coords) and only update the ligand coordinates.
    '''
    denoising_fn_from_module = functools.partial(
        self.diffusion_module, # Instance of diffusion_head.DiffusionHead
        batch=batch, # Pass the standard Batch object here
        embeddings=embeddings,
        use_conditioning=True,
    )

    num_samples = sample_config.num_samples

    # 1. Prepare masks and reference receptor coordinates
    if 'is_receptor' in batch_dict and batch_dict['is_receptor'] is not None:
        receptor_token_mask = jnp.asarray(batch_dict['is_receptor'], dtype=jnp.bool_)
    else:
        raise ValueError("No receptor mask found in batch_dict")
    
    original_complex_coords = reference_atom_coords
    atom_mask_full = reference_atom_mask

    # Tile for num_samples
    original_complex_coords_s = jnp.tile(original_complex_coords[None, ...], (num_samples, 1, 1, 1))
    atom_mask_full_s = jnp.tile(atom_mask_full[None, ...], (num_samples, 1, 1)) 

    receptor_mask_for_broadcast_s = receptor_token_mask[None, :, None, None]

    # 2. Define the per-sample denoising step for hk.scan (via hk.vmap)
    def per_sample_step_for_scan(iter_key, positions_xt_prev_step, noise_level_t_prev_step_scalar, current_noise_level_scalar_for_step):
        key_aug, current_key_noise, iter_key_next_iter = jax.random.split(iter_key, 3)

        # Augment the whole complex
        positions_xt_prev_step_augmented = diffusion_head.random_augmentation(
            rng_key=key_aug,
            positions=positions_xt_prev_step, # (num_tokens, MAX_ATOMS, 3)
            mask=atom_mask_full # (num_tokens, MAX_ATOMS)
        )

        # Noise parameters for Euler step
        gamma = sample_config.gamma_0 * (noise_level_t_prev_step_scalar > sample_config.gamma_min)
        t_hat_scalar = noise_level_t_prev_step_scalar * (1 + gamma)
        t_hat = jnp.asarray(t_hat_scalar, dtype=positions_xt_prev_step_augmented.dtype)
        
        noise_scale_val = sample_config.noise_scale * jnp.sqrt(jnp.maximum(0, t_hat_scalar**2 - noise_level_t_prev_step_scalar**2))
        noise_scale = jnp.asarray(noise_scale_val, dtype=positions_xt_prev_step_augmented.dtype)

        # Add noise to augmented structure
        added_noise_val = noise_scale * jax.random.normal(
            current_key_noise,
            positions_xt_prev_step_augmented.shape,
            dtype=positions_xt_prev_step_augmented.dtype
        )
        positions_noisy_for_denoiser = positions_xt_prev_step_augmented + added_noise_val
        
        # Core denoising call (on the whole complex)
        # The denoiser function (self.diffusion_module) returns atom positions directly
        denoised_full_structure_pred_by_af3 = denoising_fn_from_module(
            positions_noisy=positions_noisy_for_denoiser,
            noise_level=t_hat # t_hat is scalar here, or broadcastable for the denoiser
        )
        
        # Conditional Restoration: Receptor from augmented input, non-receptor from AF3 prediction
        # receptor_token_mask is (num_tokens,). Need (num_tokens, 1, 1) for per-sample broadcast
        receptor_mask_expanded_local = receptor_token_mask[:, None, None]
        
        positions_denoised_x0_hat_conditional = jnp.where(
            receptor_mask_expanded_local, # (num_tokens, 1, 1)
            positions_xt_prev_step_augmented,  # Receptor from augmented input
            denoised_full_structure_pred_by_af3    # Non-receptor from AF3's prediction
        )
        
        # Euler-Maruyama update
        safe_t_hat_scalar = jnp.where(t_hat_scalar == 0, 1e-9, t_hat_scalar)
        grad_val = (positions_noisy_for_denoiser - positions_denoised_x0_hat_conditional) / safe_t_hat_scalar
        
        d_t_val = current_noise_level_scalar_for_step - t_hat_scalar
        # This is the standard Euler update for the whole complex
        euler_updated_positions = positions_noisy_for_denoiser + sample_config.step_scale * d_t_val * grad_val

        positions_xt_next_step_final = jnp.where(
            receptor_mask_expanded_local,       # (num_tokens, 1, 1)
            positions_xt_prev_step_augmented,   # Receptor from augmented input (this step's reference)
            euler_updated_positions             # Non-receptor from the Euler update
        )

        return (iter_key_next_iter, positions_xt_next_step_final, current_noise_level_scalar_for_step)

    # Wrap per_sample_step_for_scan with hk.vmap for handling multiple samples
    # Inputs to vmapped function: iter_key_s, positions_xt_prev_step_s, noise_level_t_prev_step_s, current_noise_level_scalar
    vmapped_step_fn = hk.vmap(
        per_sample_step_for_scan,
        in_axes=(0, 0, 0, None), # key, pos, prev_noise_lvl are per sample; current_noise_lvl is shared
        out_axes=0,              # All outputs are per sample
        split_rng=False          # RNG handled manually inside
    )

    # Define the function for hk.scan
    def scan_loop_fn(carry_s, current_noise_level_scalar_for_all_samples):
        iter_key_s, positions_xt_prev_step_s, noise_level_t_prev_step_s = carry_s
        new_carry_s_tuple = vmapped_step_fn(
            iter_key_s, 
            positions_xt_prev_step_s, 
            noise_level_t_prev_step_s, 
            current_noise_level_scalar_for_all_samples
        )
        return new_carry_s_tuple, new_carry_s_tuple[1] # Return (new_carry, next_positions_s)

    # 3. Initialize noise schedule
    noise_dtype = original_complex_coords.dtype
    noise_levels_schedule = diffusion_head.noise_schedule(
        jnp.linspace(0, 1, sample_config.steps + 1)
    ).astype(noise_dtype)

    # 4. Initial positions (x_T)
    key_master_init = hk.next_rng_key()
    key_init_na_s, key_loop_init_keys_s = jax.random.split(key_master_init)
    
    if num_samples > 1:
        keys_init_na_per_sample = jax.random.split(key_init_na_s, num_samples)
    else: # Reshape to (1, key_shape)
        keys_init_na_per_sample = key_init_na_s[None, ...] 

    # Non-receptor noise generation around center_pos with specified spread
    non_receptor_noise_spread_val = 7  
    
    # center_pos is (3,). Tile for samples: (num_samples, 3)
    center_pos_s = jnp.tile(center_pos[None, :], (num_samples, 1))
    # Expand center_pos_s to (num_samples, 1, 1, 3) for broadcasting with atom-level noise
    center_pos_broadcastable_s = center_pos_s[:, None, None, :] 

    # Generate raw N(0,1) noise for all atoms in all tokens
    raw_random_displacements_s = jax.vmap(
        lambda k: jax.random.normal(k, original_complex_coords.shape, dtype=noise_dtype)
    )(keys_init_na_per_sample) # (num_samples, num_tokens, MAX_ATOMS, 3)
    
    # Scale this noise by the spread factor: N(0, non_receptor_noise_spread_val)
    scaled_random_displacements_s = raw_random_displacements_s * non_receptor_noise_spread_val 

    # Create non-receptor coordinate component: center_pos + scaled_random_displacements
    # This is applied to all tokens; jnp.where will select it for non-receptor tokens only.
    na_coords_component_s = center_pos_broadcastable_s + scaled_random_displacements_s

    # Corrected initialization of x_T:
    # Receptor part starts at its reference coordinates (unscaled by noise_levels_schedule[0])
    # Non-receptor part starts as (center + spread*noise) scaled by noise_levels_schedule[0]
    receptor_xT_s = original_complex_coords_s  # Receptor coordinates are NOT initially multiplied by sigma_max
    non_receptor_xT_s = na_coords_component_s * noise_levels_schedule[0] # Non-receptor coords (center + spread*noise) ARE scaled by sigma_max

    positions_xT_s = jnp.where(
        receptor_mask_for_broadcast_s, # (num_samples, num_tokens, 1, 1)
        receptor_xT_s,
        non_receptor_xT_s
    )
    
    # Prepare initial carry for hk.scan
    if num_samples > 1:
        scan_init_keys_for_loop_s = jax.random.split(key_loop_init_keys_s, num_samples)
    else: # Reshape to (1, key_shape)
        scan_init_keys_for_loop_s = key_loop_init_keys_s[None, ...]

    init_carry_state_s = (
        scan_init_keys_for_loop_s,                 # (num_samples, key_shape)
        positions_xT_s,                            # (num_samples, num_tokens, MAX_ATOMS, 3)
        jnp.tile(noise_levels_schedule[0], (num_samples,)), # (num_samples,) prev_noise_level_s
    )
    
    # 5. Run the scan loop
    final_carry_state_s, _ = hk.scan(
        scan_loop_fn,
        init_carry_state_s,
        noise_levels_schedule[1:], # xs for scan: t_1, t_2, ..., t_S (scalars)
        unroll=4 # As in AF3, changed from sample_config.unroll_scan_loops
    )
    
    _, final_positions_out_s, _ = final_carry_state_s # final_positions_out_s is x_0_hat_s

    # 6. Mask final results and return
    final_positions_out_masked_s = final_positions_out_s * atom_mask_full_s[..., None]

    return {'atom_positions': final_positions_out_masked_s, 'mask': atom_mask_full_s}
