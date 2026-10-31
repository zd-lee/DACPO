"""Solve DACPO's global binary orientation problem with OR-Tools CP-SAT."""

import numpy as np
from ortools.sat.python import cp_model


COST_SCALE = 1_000_000


def validate_costs(same, opposite):
    same = np.asarray(same, dtype=np.float64)
    opposite = np.asarray(opposite, dtype=np.float64)
    if same.ndim != 2 or same.shape[0] != same.shape[1]:
        raise ValueError("Expected square cost matrices")
    if same.shape != opposite.shape:
        raise ValueError("Cost matrix shapes differ")
    if not np.isfinite(same).all() or not np.isfinite(opposite).all():
        raise ValueError("Costs must be finite")
    if np.any(np.diag(same)) or np.any(np.diag(opposite)):
        raise ValueError("Self-edge costs must be zero")
    return same, opposite


def compute_cost(labels, same, opposite):
    equal = labels[:, None] == labels[None, :]
    return float(np.where(equal, same, opposite).sum())


def _solve_cp_sat(same_pair, opposite_pair, time_limit, workers):
    node_count = len(same_pair)
    if node_count == 0:
        return np.zeros(0, dtype=np.int32), {
            "status": "empty",
            "wall_time": 0.0,
            "integer_best_bound": 0.0,
        }

    model = cp_model.CpModel()
    flips = [model.new_bool_var(f"flip_{idx}") for idx in range(node_count)]
    # Flipping every label leaves the objective unchanged.
    model.add(flips[0] == 0)

    objective_terms = []
    rows, cols = np.triu_indices(node_count, k=1)
    for edge_idx, (src, dst) in enumerate(zip(rows, cols)):
        delta = opposite_pair[src, dst] - same_pair[src, dst]
        coefficient = int(round(delta * COST_SCALE))
        if coefficient == 0:
            continue
        different = model.new_bool_var(f"different_{edge_idx}")
        model.add(flips[src] + flips[dst] == 1).only_enforce_if(different)
        model.add(flips[src] == flips[dst]).only_enforce_if(different.Not())
        objective_terms.append(coefficient * different)

    model.minimize(sum(objective_terms))
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit
    solver.parameters.num_search_workers = workers
    status = solver.solve(model)
    if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        raise RuntimeError(f"CP-SAT did not produce a solution: {solver.status_name(status)}")

    labels = np.array([solver.value(var) for var in flips], dtype=np.int32)
    return labels, {
        "status": solver.status_name(status),
        "wall_time": solver.wall_time,
        "integer_best_bound": solver.best_objective_bound,
    }


def solve_flip(same, opposite, time_limit=300.0, workers=4, allow_feasible=False):
    same, opposite = validate_costs(same, opposite)
    if not np.isfinite(time_limit) or time_limit <= 0:
        raise ValueError("time_limit must be positive and finite")
    if isinstance(workers, bool) or not isinstance(workers, int) or not 1 <= workers <= 32:
        raise ValueError("workers must be an integer in [1, 32]")

    # DACPO stores both directed entries. Aggregate them once per undirected edge.
    same_pair = same + same.T
    opposite_pair = opposite + opposite.T
    labels, solver_stats = _solve_cp_sat(
        same_pair, opposite_pair, time_limit=time_limit, workers=workers
    )
    if not allow_feasible and solver_stats["status"] not in {"OPTIMAL", "empty"}:
        raise RuntimeError(f"CP-SAT has not proved optimality: {solver_stats}")
    if labels.shape != (len(same),) or not np.isin(labels, [0, 1]).all():
        raise RuntimeError("CP-SAT returned invalid binary labels")

    rows, cols = np.triu_indices(len(same), k=1)
    deltas = opposite_pair[rows, cols] - same_pair[rows, cols]
    active_edges = int(np.count_nonzero(deltas))
    rounding_half_width = active_edges / (2 * COST_SCALE)
    constant_cost = float(same_pair[rows, cols].sum())
    cost = compute_cost(labels, same, opposite)
    stats = {
        "method": "cp-sat",
        "status": solver_stats["status"],
        "nodes": len(same),
        "edges": active_edges,
        "objective": cost,
        "dacpo_cost": cost,
        "initial_cost": float(same.sum()),
        "wall_time": solver_stats["wall_time"],
        "cost_scale": COST_SCALE,
        "cost_quantization_bound": active_edges / COST_SCALE,
        "cost_lower_bound": (
            constant_cost
            + solver_stats["integer_best_bound"] / COST_SCALE
            - rounding_half_width
        ),
    }
    return labels, stats
