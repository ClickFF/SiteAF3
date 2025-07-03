import numpy as np


def calculate_ilDDT(
    pred_pos,  # (B, L, 3) or (L, 3)
    true_pos,  # (B, L, 3) or (L, 3)
    interface_mask,  # (B, L) or (L,)
    per_atom=False,  # if true, return lDDT for each atom (retain L axis)
    true_dist_cutoff=15.0,
    dist_diff_cutoffs=[0.5, 1.0, 2.0, 4.0],
):
    # Calculate predicted and true distance matrices
    def pos_to_dist(pos):
        return ((pos[..., None, :, :] - pos[..., None, :]) ** 2).sum(-1) ** 0.5

    pred_dist, true_dist = pos_to_dist(pred_pos), pos_to_dist(true_pos)

    # Get mask for distances to be scored (remove self distances)
    L = pred_pos.shape[-2]
    dist_mask = (true_dist < true_dist_cutoff) * (1.0 - np.eye(L))

    # Calculate contribution of each distance to lDDT
    dist_diff = np.abs(pred_dist - true_dist)
    score = (dist_diff[..., None] < np.array(dist_diff_cutoffs)).mean(-1)

    # Return calculated lDDT
    axis = (-1,) if per_atom else (-2, -1)
    return (score * dist_mask).sum(axis) / dist_mask.sum(axis)