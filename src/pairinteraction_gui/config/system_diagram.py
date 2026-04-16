# SPDX-FileCopyrightText: 2026 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

from math import cos, hypot, radians, sin
from typing import TYPE_CHECKING, Any

from PySide6.QtCore import QPointF, Qt
from PySide6.QtGui import QColor, QFont, QPainter, QPen, QPolygonF
from PySide6.QtWidgets import QSizePolicy, QWidget

if TYPE_CHECKING:
    from PySide6.QtGui import QPaintEvent

    from pairinteraction_gui.config.system_config import SystemConfig


class _Colors:
    AXIS = QColor("#8a91a0")
    ATOM_FILL = QColor("#1a1a1a")
    ATOM_LABEL = QColor("#ffffff")
    ION_POS = QColor("#cc2222")
    ION_NEG = QColor("#2266cc")
    ION_NEUTRAL = QColor("#8a91a0")
    EFIELD = QColor("#e05050")
    BFIELD = QColor("#2266cc")
    TEXT = QColor("#111828")
    BG = QColor("#f9fbfe")


class SystemDiagram(QWidget):
    """Schematic diagram showing atom(s), fields, and ion position."""

    def __init__(self, parent: QWidget, *, two_atoms: bool) -> None:
        super().__init__(parent)
        self._two_atoms = two_atoms
        self._config: SystemConfig | None = None
        self.setMinimumHeight(120)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)

    def connectConfig(self, config: SystemConfig) -> None:
        """Connect all relevant config items to trigger repaints."""
        for item in [
            config.Ex,
            config.Ey,
            config.Ez,
            config.Bx,
            config.By,
            config.Bz,
            config.ion_distance,
            config.ion_angle,
            config.ion_charge,
        ]:
            item.connectAll(self.update)
        if self._two_atoms:
            config.distance.connectAll(self.update)  # type: ignore[attr-defined]
            config.angle.connectAll(self.update)  # type: ignore[attr-defined]
        self._config = config

    def _read_state(self, use_max: bool) -> dict[str, Any]:
        cfg = self._config
        if cfg is None:
            return {}

        def spinbox(item: Any) -> Any:
            return item.max_spinbox if use_max else item.min_spinbox

        ex = spinbox(cfg.Ex).value() if cfg.Ex.isChecked() else 0.0
        ey = spinbox(cfg.Ey).value() if cfg.Ey.isChecked() else 0.0
        ez = spinbox(cfg.Ez).value() if cfg.Ez.isChecked() else 0.0
        bx = spinbox(cfg.Bx).value() if cfg.Bx.isChecked() else 0.0
        by = spinbox(cfg.By).value() if cfg.By.isChecked() else 0.0
        bz = spinbox(cfg.Bz).value() if cfg.Bz.isChecked() else 0.0

        def _sign(v: float) -> int:
            return 1 if v > 1e-9 else (-1 if v < -1e-9 else 0)

        state: dict[str, Any] = {
            "has_efield": hypot(ex, ez) > 1e-9,
            "efield_dir": (ex, ez),
            "ey_sign": _sign(ey),
            "has_bfield": hypot(bx, bz) > 1e-9,
            "bfield_dir": (bx, bz),
            "by_sign": _sign(by),
            "show_ion": cfg.ion_distance.isChecked(),
            "ion_angle_deg": spinbox(cfg.ion_angle).value() if cfg.ion_angle.isChecked() else 0.0,
            "ion_distance_val": spinbox(cfg.ion_distance).value() if cfg.ion_distance.isChecked() else 0.0,
            "ion_charge": cfg.ion_charge.spinbox.value(),
        }

        if self._two_atoms:
            state["atom_angle_deg"] = spinbox(cfg.angle).value()  # type: ignore[attr-defined]
            dist_item = cfg.distance  # type: ignore[attr-defined]
            state["atom_distance_val"] = spinbox(dist_item).value() if dist_item.isChecked() else None

        return state

    def paintEvent(self, event: QPaintEvent) -> None:
        if self._config is None:
            return

        w, h = self.width(), self.height()
        half = w // 2

        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.fillRect(0, 0, w, h, _Colors.BG)

        painter.setPen(QPen(_Colors.AXIS, 1, Qt.PenStyle.DotLine))
        painter.drawLine(half, 4, half, h - 4)

        for use_max, x_offset, label in [(False, 0, "min"), (True, half, "max")]:
            painter.save()
            painter.setClipRect(x_offset, 0, half, h)
            self._draw_panel(painter, x_offset, half, h, self._read_state(use_max), label)
            painter.restore()

        painter.end()

    def _draw_panel(
        self,
        painter: QPainter,
        x_offset: int,
        panel_w: int,
        panel_h: int,
        state: dict[str, Any],
        label: str,
    ) -> None:
        cy = panel_h / 2.0
        origin_px = QPointF(x_offset + panel_w * 0.22, cy)
        atoms_cx = x_offset + panel_w * 0.68
        scale = min(panel_w * 0.5, panel_h) * 0.30

        font = QFont()
        font.setPointSize(7)
        painter.setFont(font)
        painter.setPen(QPen(_Colors.AXIS))
        painter.drawText(
            int(x_offset),
            2,
            panel_w,
            14,
            Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignTop,
            label,
        )

        if self._two_atoms:
            angle_rad = radians(state.get("atom_angle_deg", 0.0))
            dx, dz = sin(angle_rad) * scale, cos(angle_rad) * scale
            atom1_px = QPointF(atoms_cx - dx / 2, cy + dz / 2)
            atom2_px: QPointF | None = QPointF(atoms_cx + dx / 2, cy - dz / 2)
        else:
            atom1_px = QPointF(atoms_cx, cy)
            atom2_px = None

        ion_px = None
        if state.get("show_ion"):
            ion_rad = radians(state.get("ion_angle_deg", 0.0))
            atom_dist = state.get("atom_distance_val")  # None means infinite
            if not self._two_atoms:
                ion_scale = scale
            elif self._two_atoms and atom_dist:
                ion_dist = state.get("ion_distance_val", 0.0)
                ion_scale = (ion_dist / atom_dist) * scale
            else:
                ion_scale = 0.1 * scale
            ion_px = QPointF(
                atom1_px.x() + sin(ion_rad) * ion_scale,
                atom1_px.y() - cos(ion_rad) * ion_scale,
            )

        self._draw_axes(painter, origin_px)

        if atom2_px is not None:
            self._draw_atom(painter, atom2_px, "2")
        self._draw_atom(painter, atom1_px, "1" if self._two_atoms else "")

        if ion_px is not None:
            self._draw_ion(painter, ion_px, state.get("ion_charge", 1.0))

        if state.get("has_efield"):
            ex, ez = state["efield_dir"]
            mag = hypot(ex, ez)
            self._draw_arrow(painter, origin_px, (ex / mag, ez / mag), scale, _Colors.EFIELD, "E")

        if state.get("has_bfield"):
            bx, bz = state["bfield_dir"]
            mag = hypot(bx, bz)
            self._draw_arrow(painter, origin_px, (bx / mag, bz / mag), scale, _Colors.BFIELD, "B")

        # Out-of-plane indicators stacked below the axis origin
        oop_y = origin_px.y() + 18
        for sign_key, color, field_label in [
            ("ey_sign", _Colors.EFIELD, "E"),
            ("by_sign", _Colors.BFIELD, "B"),
        ]:
            if state.get(sign_key, 0) != 0:
                self._draw_oop_symbol(painter, QPointF(origin_px.x(), oop_y), state[sign_key], color, field_label)
                oop_y += 16

    def _draw_axes(self, painter: QPainter, origin: QPointF) -> None:
        axis_len = 20.0
        painter.setPen(QPen(_Colors.AXIS, 1.2, Qt.PenStyle.DashLine))
        x_end = QPointF(origin.x() + axis_len, origin.y())
        z_end = QPointF(origin.x(), origin.y() - axis_len)
        painter.drawLine(origin, x_end)
        painter.drawLine(origin, z_end)

        font = QFont()
        font.setPointSize(7)
        painter.setFont(font)
        painter.setPen(QPen(_Colors.AXIS))
        painter.drawText(x_end + QPointF(2, 4), "x")
        painter.drawText(z_end + QPointF(2, -2), "z")

    def _draw_atom(self, painter: QPainter, center: QPointF, label: str) -> None:
        r = 8.0
        painter.setBrush(_Colors.ATOM_FILL)
        painter.setPen(Qt.PenStyle.NoPen)
        painter.drawEllipse(center, r, r)

        if label:
            font = QFont()
            font.setPointSize(7)
            font.setBold(True)
            painter.setFont(font)
            painter.setPen(QPen(_Colors.ATOM_LABEL))
            painter.drawText(
                int(center.x() - r),
                int(center.y() - r),
                int(r * 2),
                int(r * 2),
                Qt.AlignmentFlag.AlignCenter,
                label,
            )

    def _draw_ion(self, painter: QPainter, center: QPointF, charge: float) -> None:
        r = 5.0
        if charge > 1e-9:
            color = _Colors.ION_POS
            sign = "+"
        elif charge < -1e-9:
            color = _Colors.ION_NEG
            sign = "\u2212"  # minus sign
        else:
            color = _Colors.ION_NEUTRAL
            sign = "0"

        painter.setBrush(color)
        painter.setPen(Qt.PenStyle.NoPen)
        painter.drawEllipse(center, r, r)

        inner = r * 0.55
        painter.setPen(QPen(QColor("#ffffff"), 1.5))
        painter.drawLine(QPointF(center.x() - inner, center.y()), QPointF(center.x() + inner, center.y()))
        if sign == "+":
            painter.drawLine(QPointF(center.x(), center.y() - inner), QPointF(center.x(), center.y() + inner))

    def _draw_arrow(
        self,
        painter: QPainter,
        origin: QPointF,
        direction: tuple[float, float],
        scale: float,
        color: QColor,
        label: str,
    ) -> None:
        arrow_len = min(scale * 0.5, 40.0)
        px, pz = direction  # physics x, z; pixel y is inverted
        tip = QPointF(origin.x() + px * arrow_len, origin.y() - pz * arrow_len)

        painter.setPen(QPen(color, 1.8))
        painter.drawLine(origin, tip)

        painter.setBrush(color)
        painter.setPen(Qt.PenStyle.NoPen)
        painter.drawPolygon(self._make_arrowhead(tip, (px, -pz)))

        font = QFont()
        font.setPointSize(7)
        font.setItalic(True)
        painter.setFont(font)
        painter.setPen(QPen(color))
        painter.drawText(tip + QPointF(px * 5 + 2, -pz * 5 + 4), label)

    @staticmethod
    def _draw_oop_symbol(
        painter: QPainter,
        center: QPointF,
        sign: int,
        color: QColor,
        label: str,
    ) -> None:
        """Draw ⊙ (out-of-page, sign=+1) or ⊗ (into-page, sign=-1) with a field label."""
        r = 5.0
        painter.setBrush(Qt.BrushStyle.NoBrush)
        painter.setPen(QPen(color, 1.2))
        painter.drawEllipse(center, r, r)

        if sign > 0:
            painter.setBrush(color)
            painter.setPen(Qt.PenStyle.NoPen)
            painter.drawEllipse(center, 1.5, 1.5)
        else:
            d = r * 0.6
            painter.setPen(QPen(color, 1.2))
            painter.drawLine(QPointF(center.x() - d, center.y() - d), QPointF(center.x() + d, center.y() + d))
            painter.drawLine(QPointF(center.x() + d, center.y() - d), QPointF(center.x() - d, center.y() + d))

        font = QFont()
        font.setPointSize(7)
        font.setItalic(True)
        painter.setFont(font)
        painter.setPen(QPen(color))
        painter.drawText(center + QPointF(r + 2, 4), label)

    @staticmethod
    def _make_arrowhead(tip: QPointF, direction: tuple[float, float]) -> QPolygonF:
        ux, uy = direction
        px, py = -uy, ux  # perpendicular
        head_len, head_w = 8.0, 4.0
        base = QPointF(tip.x() - head_len * ux, tip.y() - head_len * uy)
        return QPolygonF(
            [
                tip,
                QPointF(base.x() + head_w * px, base.y() + head_w * py),
                QPointF(base.x() - head_w * px, base.y() - head_w * py),
            ]
        )
