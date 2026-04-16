# SPDX-FileCopyrightText: 2026 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

from math import cos, hypot, radians, sin
from typing import TYPE_CHECKING

from PySide6.QtCore import QPointF, Qt
from PySide6.QtGui import QColor, QFont, QPainter, QPen, QPolygonF
from PySide6.QtWidgets import QSizePolicy, QWidget

if TYPE_CHECKING:
    from pairinteraction_gui.config.system_config import SystemConfig


class _Colors:
    AXIS = QColor("#8a91a0")
    ATOM_FILL = QColor("#cad9f7")
    ATOM_EDGE = QColor("#28354e")
    ION_FILL = QColor("#f9fbfe")
    ION_EDGE = QColor("#28354e")
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
        for item in [config.Ex, config.Ey, config.Ez, config.Bx, config.By, config.Bz,
                     config.ion_distance, config.ion_angle]:
            item.connectAll(self.update)
        if self._two_atoms:
            config.distance.connectAll(self.update)  # type: ignore[attr-defined]
            config.angle.connectAll(self.update)  # type: ignore[attr-defined]
        self._config = config

    def _read_state(self, use_max: bool) -> dict:
        cfg = self._config
        if cfg is None:
            return {}

        spinbox = lambda item: item.max_spinbox if use_max else item.min_spinbox  # noqa: E731

        ex = spinbox(cfg.Ex).value() if cfg.Ex.isChecked() else 0.0
        ez = spinbox(cfg.Ez).value() if cfg.Ez.isChecked() else 0.0
        bx = spinbox(cfg.Bx).value() if cfg.Bx.isChecked() else 0.0
        bz = spinbox(cfg.Bz).value() if cfg.Bz.isChecked() else 0.0
        # Ey/By out-of-plane: read but not drawn

        has_efield = hypot(ex, ez) > 1e-9
        has_bfield = hypot(bx, bz) > 1e-9

        show_ion = cfg.ion_distance.isChecked()
        ion_angle_deg = spinbox(cfg.ion_angle).value() if cfg.ion_angle.isChecked() else 0.0

        state: dict = {
            "has_efield": has_efield,
            "efield_dir": (ex, ez),
            "has_bfield": has_bfield,
            "bfield_dir": (bx, bz),
            "show_ion": show_ion,
            "ion_angle_deg": ion_angle_deg,
        }

        if self._two_atoms:
            state["atom_angle_deg"] = spinbox(cfg.angle).value()  # type: ignore[attr-defined]

        return state

    def paintEvent(self, event) -> None:  # noqa: ANN001
        w, h = self.width(), self.height()
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)

        painter.fillRect(0, 0, w, h, _Colors.BG)

        half = w // 2

        # Draw dividing line
        painter.setPen(QPen(_Colors.AXIS, 1, Qt.PenStyle.DotLine))
        painter.drawLine(half, 4, half, h - 4)

        # Left panel: min step
        painter.save()
        painter.setClipRect(0, 0, half, h)
        self._draw_panel(painter, 0, half, h, self._read_state(use_max=False), label="min")
        painter.restore()

        # Right panel: max step
        painter.save()
        painter.setClipRect(half, 0, half, h)
        self._draw_panel(painter, half, half, h, self._read_state(use_max=True), label="max")
        painter.restore()

        painter.end()

    def _draw_panel(
        self,
        painter: QPainter,
        x_offset: int,
        panel_w: int,
        panel_h: int,
        state: dict,
        label: str,
    ) -> None:
        cy = panel_h / 2.0

        origin_px = QPointF(x_offset + panel_w * 0.22, cy)
        atoms_cx = x_offset + panel_w * 0.68
        scale = min(panel_w * 0.5, panel_h) * 0.30

        # Panel label ("min" / "max") at top-center
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

        if self._two_atoms and state:
            atom_angle_deg = state.get("atom_angle_deg", 0.0)
            angle_rad = radians(atom_angle_deg)
            dx = sin(angle_rad) * scale
            dz = cos(angle_rad) * scale

            atom1_px = QPointF(atoms_cx - dx / 2, cy + dz / 2)
            atom2_px = QPointF(atoms_cx + dx / 2, cy - dz / 2)
        else:
            atom1_px = QPointF(atoms_cx, cy)
            atom2_px = None

        # Ion position (relative to atom1)
        ion_px = None
        if state.get("show_ion"):
            ion_angle_rad = radians(state.get("ion_angle_deg", 0.0))
            ion_dx = sin(ion_angle_rad) * scale * 0.7
            ion_dz = cos(ion_angle_rad) * scale * 0.7
            ion_px = QPointF(atom1_px.x() + ion_dx, atom1_px.y() - ion_dz)

        self._draw_axes(painter, origin_px)

        if atom2_px is not None:
            self._draw_atom(painter, atom2_px, "2")
        self._draw_atom(painter, atom1_px, "1" if self._two_atoms else "")

        if ion_px is not None:
            self._draw_ion(painter, ion_px)

        if state.get("has_efield"):
            ex, ez = state["efield_dir"]
            mag = hypot(ex, ez)
            self._draw_arrow(painter, origin_px, (ex / mag, ez / mag), scale, _Colors.EFIELD, "E")

        if state.get("has_bfield"):
            bx, bz = state["bfield_dir"]
            mag = hypot(bx, bz)
            self._draw_arrow(painter, origin_px, (bx / mag, bz / mag), scale, _Colors.BFIELD, "B")

    def _draw_axes(self, painter: QPainter, origin: QPointF) -> None:
        axis_len = 20.0
        pen = QPen(_Colors.AXIS, 1.2, Qt.PenStyle.DashLine)
        painter.setPen(pen)

        # X axis (right)
        x_end = QPointF(origin.x() + axis_len, origin.y())
        painter.drawLine(origin, x_end)

        # Z axis (up)
        z_end = QPointF(origin.x(), origin.y() - axis_len)
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
        painter.setPen(QPen(_Colors.ATOM_EDGE, 1.5))
        painter.drawEllipse(center, r, r)

        if label:
            font = QFont()
            font.setPointSize(7)
            painter.setFont(font)
            painter.setPen(QPen(_Colors.TEXT))
            painter.drawText(
                int(center.x() - r),
                int(center.y() - r),
                int(r * 2),
                int(r * 2),
                Qt.AlignmentFlag.AlignCenter,
                label,
            )

    def _draw_ion(self, painter: QPainter, center: QPointF) -> None:
        r = 4.0
        pen = QPen(_Colors.ION_EDGE, 1.5)
        painter.setPen(pen)
        painter.drawLine(
            QPointF(center.x() - r, center.y()),
            QPointF(center.x() + r, center.y()),
        )
        painter.drawLine(
            QPointF(center.x(), center.y() - r),
            QPointF(center.x(), center.y() + r),
        )

        font = QFont()
        font.setPointSize(7)
        painter.setFont(font)
        painter.setPen(QPen(_Colors.TEXT))
        painter.drawText(center + QPointF(r + 2, 4), "ion")

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
        dx_phys, dz_phys = direction
        # pixel direction: x right, y down → z up means py = -dz_phys
        pdx = dx_phys
        pdy = -dz_phys

        tip = QPointF(origin.x() + pdx * arrow_len, origin.y() + pdy * arrow_len)

        pen = QPen(color, 1.8)
        painter.setPen(pen)
        painter.drawLine(origin, tip)

        # arrowhead
        arrowhead = self._make_arrowhead(tip, (pdx, pdy))
        painter.setBrush(color)
        painter.setPen(Qt.PenStyle.NoPen)
        painter.drawPolygon(arrowhead)

        # label
        font = QFont()
        font.setPointSize(7)
        font.setItalic(True)
        painter.setFont(font)
        painter.setPen(QPen(color))
        painter.drawText(tip + QPointF(pdx * 5 + 2, pdy * 5 + 4), label)

    @staticmethod
    def _make_arrowhead(tip: QPointF, direction: tuple[float, float]) -> QPolygonF:
        ux, uy = direction
        # perpendicular
        px, py = -uy, ux
        L = 8.0
        W = 4.0
        base = QPointF(tip.x() - L * ux, tip.y() - L * uy)
        left = QPointF(base.x() + W * px, base.y() + W * py)
        right = QPointF(base.x() - W * px, base.y() - W * py)
        return QPolygonF([tip, left, right])
