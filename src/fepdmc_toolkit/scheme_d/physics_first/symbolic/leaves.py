"""S02 unique public IR writer for symbolic leaves. No numeric identity."""

from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np

from fepdmc_toolkit.scheme_d.ir import Binding, Layout

_NUMERIC_NAME = re.compile(r"(?:\d\.\d|[eE][+-]\d|\.16e)")


def assert_symbolic_name(name: str) -> str:
    if _NUMERIC_NAME.search(name) or "." in name:
        raise ValueError(f"numeric identity forbidden in leaf name: {name}")
    return name


@dataclass(frozen=True)
class IntegerMomentumForm:
    k_coeff: int
    q_coeffs: tuple[tuple[int, int], ...]

    @classmethod
    def electron_open(cls, open_slots: tuple[int, ...]) -> IntegerMomentumForm:
        slots = tuple(sorted(set(int(s) for s in open_slots)))
        return cls(1, tuple((s, -1) for s in slots))

    def evaluate(self, k: float, q: np.ndarray) -> float:
        acc = self.k_coeff * float(k)
        for slot, coeff in self.q_coeffs:
            acc += coeff * float(q[slot])
        return acc

    def tag(self) -> str:
        if not self.q_coeffs:
            return f"k{self.k_coeff}"
        qpart = "_".join(f"s{s}c{c}" for s, c in self.q_coeffs)
        return f"k{self.k_coeff}_{qpart}"


@dataclass(frozen=True)
class ElectronProp:
    provider: str
    open_q_slots: tuple[int, ...]
    tau_left_idx: int
    tau_right_idx: int
    electron_interface: str
    momentum: IntegerMomentumForm

    def input_name(self) -> str:
        slots = "none" if not self.open_q_slots else "_".join(str(s) for s in self.open_q_slots)
        return assert_symbolic_name(
            f"E_p{self.provider}_Q{slots}_T{self.tau_left_idx}_{self.tau_right_idx}_I{self.electron_interface}"
        )

    def runtime_dependencies(self) -> tuple[str, ...]:
        deps = ["k", "t", f"tau[{self.tau_left_idx}]", f"tau[{self.tau_right_idx}]"]
        deps.extend(f"q[{s}]" for s in self.open_q_slots)
        if self.electron_interface == "matrix2":
            deps.extend(("delta", "gap"))
        return tuple(deps)

    def layout_shape(self) -> tuple[str, tuple[int, int]]:
        if self.electron_interface == "scalar":
            return Layout.SCALAR.value, (1, 1)
        if self.electron_interface == "matrix2":
            return Layout.DENSE.value, (2, 2)
        raise ValueError(self.electron_interface)


@dataclass(frozen=True)
class PhononProp:
    provider: str
    q_slot: int
    mode_slot: int
    tau_open_idx: int
    tau_close_idx: int
    layout_bands: int = 1

    def input_name(self) -> str:
        return assert_symbolic_name(
            f"D_p{self.provider}_Q{self.q_slot}_M{self.mode_slot}_T{self.tau_open_idx}_{self.tau_close_idx}_B{self.layout_bands}"
        )

    def runtime_dependencies(self) -> tuple[str, ...]:
        return (
            f"q[{self.q_slot}]",
            f"tau[{self.tau_open_idx}]",
            f"tau[{self.tau_close_idx}]",
            "omega",
        )

    def layout_shape(self) -> tuple[str, tuple[int, int]]:
        if self.layout_bands == 1:
            return Layout.SCALAR.value, (1, 1)
        b = self.layout_bands
        return Layout.DIAGONAL.value, (b, b)


@dataclass(frozen=True)
class VertexProp:
    provider: str
    direction: str
    incoming_momentum: IntegerMomentumForm
    q_slot: int
    mode_slot: int
    band_interface: str

    def input_name(self) -> str:
        return assert_symbolic_name(
            f"V_p{self.provider}_D{self.direction}_Q{self.q_slot}_M{self.mode_slot}_I{self.band_interface}_F{self.incoming_momentum.tag()}"
        )

    def runtime_dependencies(self) -> tuple[str, ...]:
        deps = ["g", "L", f"q[{self.q_slot}]"]
        if self.incoming_momentum.q_coeffs:
            deps.extend(f"q[{s}]" for s, _ in self.incoming_momentum.q_coeffs)
            deps.append("k")
        return tuple(deps)

    def layout_shape(self) -> tuple[str, tuple[int, int]]:
        if self.band_interface == "scalar":
            return Layout.SCALAR.value, (1, 1)
        if self.band_interface == "matrix2":
            return Layout.DIAGONAL.value, (2, 2)
        raise ValueError(self.band_interface)


@dataclass(frozen=True)
class PrefactorProp:
    n: int
    normalization_spec: str

    def input_name(self) -> str:
        return assert_symbolic_name(f"PF_n{self.n}_N{self.normalization_spec}")

    def runtime_dependencies(self) -> tuple[str, ...]:
        return ("g", "L")

    def layout_shape(self) -> tuple[str, tuple[int, int]]:
        return Layout.SCALAR.value, (1, 1)


@dataclass(frozen=True)
class BoundaryProp:
    boundary_type: str
    external_index_slots: tuple[int, ...]
    provider: str

    def input_name(self) -> str:
        slots = "_".join(str(s) for s in self.external_index_slots) or "none"
        return assert_symbolic_name(f"B_p{self.provider}_T{self.boundary_type}_X{slots}")

    def runtime_dependencies(self) -> tuple[str, ...]:
        return tuple(f"ext[{s}]" for s in self.external_index_slots)

    def layout_shape(self) -> tuple[str, tuple[int, int]]:
        return Layout.SCALAR.value, (1, 1)


LeafDesc = ElectronProp | PhononProp | VertexProp | PrefactorProp | BoundaryProp


def electron_prop(
    *,
    provider: str,
    open_q_slots: tuple[int, ...],
    tau_left_idx: int,
    tau_right_idx: int,
    electron_interface: str,
) -> ElectronProp:
    slots = tuple(sorted(int(s) for s in open_q_slots))
    if tau_left_idx == tau_right_idx:
        raise ValueError("electron time endpoints must be distinct indices")
    return ElectronProp(
        provider=provider,
        open_q_slots=slots,
        tau_left_idx=int(tau_left_idx),
        tau_right_idx=int(tau_right_idx),
        electron_interface=electron_interface,
        momentum=IntegerMomentumForm.electron_open(slots),
    )


def phonon_prop(
    *,
    provider: str,
    q_slot: int,
    mode_slot: int,
    tau_open_idx: int,
    tau_close_idx: int,
    layout_bands: int = 1,
) -> PhononProp:
    if tau_open_idx == tau_close_idx:
        raise ValueError("phonon time endpoints must be distinct indices")
    return PhononProp(
        provider=provider,
        q_slot=int(q_slot),
        mode_slot=int(mode_slot),
        tau_open_idx=int(tau_open_idx),
        tau_close_idx=int(tau_close_idx),
        layout_bands=int(layout_bands),
    )


def vertex_prop(
    *,
    provider: str,
    direction: str,
    incoming_open_slots: tuple[int, ...],
    q_slot: int,
    mode_slot: int,
    band_interface: str,
) -> VertexProp:
    if direction not in {"emit", "absorb"}:
        raise ValueError(direction)
    return VertexProp(
        provider=provider,
        direction=direction,
        incoming_momentum=IntegerMomentumForm.electron_open(incoming_open_slots),
        q_slot=int(q_slot),
        mode_slot=int(mode_slot),
        band_interface=band_interface,
    )


def prefactor_prop(*, n: int, normalization_spec: str) -> PrefactorProp:
    return PrefactorProp(n=int(n), normalization_spec=normalization_spec)


class LeafTable:
    """Stable ABI order. Unique writer for names used as Graph input var_name."""

    def __init__(self) -> None:
        self.order: list[str] = []
        self.descriptors: dict[str, LeafDesc] = {}

    def register(self, desc: LeafDesc) -> str:
        name = desc.input_name()
        existing = self.descriptors.get(name)
        if existing is None:
            self.descriptors[name] = desc
            self.order.append(name)
            return name
        if existing != desc:
            raise ValueError(f"leaf name collision {name}")
        return name

    def __len__(self) -> int:
        return len(self.order)

    def rows(self) -> list[dict]:
        out = []
        for name in self.order:
            d = self.descriptors[name]
            layout, shape = d.layout_shape()
            out.append(
                {
                    "name": name,
                    "kind": type(d).__name__,
                    "layout": layout,
                    "shape": "x".join(str(s) for s in shape),
                    "runtime_dependencies": ";".join(d.runtime_dependencies()),
                    "forbidden_identity": "numeric momenta/time or float strings",
                }
            )
        return out


def _scalar_arr(z: complex) -> np.ndarray:
    return np.array([[complex(z)]], dtype=np.complex128)


def eval_leaf(desc: LeafDesc, binding) -> np.ndarray:
    """Numeric value of one symbolic leaf. No graph rebuild."""
    q = np.asarray(binding.q, dtype=np.float64)
    tau = np.asarray(binding.tau, dtype=np.float64)
    if isinstance(desc, PrefactorProp):
        if desc.normalization_spec != "g_over_sqrt_L":
            raise ValueError(desc.normalization_spec)
        gv = binding.g / np.sqrt(binding.L)
        return _scalar_arr(gv ** (2 * desc.n))
    if isinstance(desc, ElectronProp):
        k_eff = desc.momentum.evaluate(binding.k, q)
        dtau = float(tau[desc.tau_right_idx] - tau[desc.tau_left_idx])
        xi = 2.0 * binding.t * (1.0 - np.cos(k_eff))
        if desc.electron_interface == "scalar":
            return _scalar_arr(np.exp(-xi * dtau))
        if desc.electron_interface == "matrix2":
            h00 = xi
            h11 = xi + binding.gap
            h01 = binding.delta
            H = np.array([[h00, h01], [h01, h11]], dtype=np.complex128)
            w, v = np.linalg.eigh(H)
            return (v * np.exp(-w * dtau)) @ v.conj().T
        raise ValueError(desc.electron_interface)
    if isinstance(desc, PhononProp):
        dtau = float(tau[desc.tau_close_idx] - tau[desc.tau_open_idx])
        dval = np.exp(-binding.omega * dtau)
        if desc.layout_bands == 1:
            return _scalar_arr(dval)
        eye = np.eye(desc.layout_bands, dtype=np.complex128)
        return eye * dval
    if isinstance(desc, VertexProp):
        gv = binding.g / np.sqrt(binding.L)
        if desc.band_interface == "scalar":
            return _scalar_arr(gv)
        if desc.band_interface == "matrix2":
            return gv * np.array([[1.0, 0.0], [0.0, -1.0]], dtype=np.complex128)
        raise ValueError(desc.band_interface)
    if isinstance(desc, BoundaryProp):
        raise ValueError("current object has no external boundary leaf")
    raise TypeError(type(desc))


def graph_input(builder, table: LeafTable, desc: LeafDesc) -> int:
    name = table.register(desc)
    layout, shape = desc.layout_shape()
    return builder.input(name, shape, layout, Binding.BOUND.value)
