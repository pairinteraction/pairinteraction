# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

import logging
import os
import signal
from typing import TYPE_CHECKING

from PySide6.QtCore import QObject, QRectF, QSocketNotifier, Qt, QTimer, Signal
from PySide6.QtGui import QColor, QFont, QPainter, QPainterPath, QPen, QPixmap
from PySide6.QtWidgets import QApplication, QSplashScreen

from pairinteraction_gui.worker import MultiThreadWorker

if TYPE_CHECKING:
    from types import FrameType

logger = logging.getLogger(__name__)


class MainSignals(QObject):
    """Signals for the application.

    We store an instance of this signal class in the Application instance, see app.py.
    So to access these signals (from anywhere in the application), you can use
    `Application.instance().signals`.
    """

    ask_download_database = Signal(str)


class Application(QApplication):
    """Add some global signals to the QApplication."""

    signals = MainSignals()

    @staticmethod
    def instance() -> Application:
        """Return the current instance of the application."""
        return QApplication.instance()  # type: ignore [return-value]

    def allow_ctrl_c(self) -> None:
        # Create a pipe to communicate between the signal handler and the Qt event loop
        pipe_r, pipe_w = os.pipe()

        def signal_handler(signal: int, frame: FrameType | None) -> None:
            os.write(pipe_w, b"x")  # Write a single byte to the pipe

        signal.signal(signal.SIGINT, signal_handler)

        def handle_signal() -> None:
            os.read(pipe_r, 1)  # Read the byte from the pipe to clear it
            logger.info("Ctrl+C detected in terminal. Shutting down gracefully...")
            self.quit()

        sn = QSocketNotifier(pipe_r, QSocketNotifier.Type.Read, parent=self)
        sn.activated.connect(handle_signal)

        # Create a timer to ensure the event loop processes events regularly
        # This makes Ctrl+C work even when the application is idle
        timer = QTimer(self)
        timer.timeout.connect(lambda: None)  # Do nothing, just wake up the event loop
        timer.start(200)

    @staticmethod
    def quit() -> None:
        """Quit the application."""
        logger.debug("Calling Application.quit().")
        MultiThreadWorker.terminate_all()
        QApplication.quit()
        logger.debug("Application.quit() done.")


class SplashScreen(QSplashScreen):
    """Startup splash screen with the PairInteraction logo."""

    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("Launching the PairInteraction GUI")
        border = 16
        logo_size = 250
        full_size = logo_size + 2 * border
        pixmap = QPixmap(full_size, full_size)
        pixmap.fill(QColor("#28354e"))  # "Dark" theme color
        painter = QPainter(pixmap)
        painter.translate(border, border)
        _draw_logo(painter, logo_size)
        painter.end()
        self.setPixmap(pixmap)


def _draw_logo(painter: QPainter, size: int) -> None:
    """Draw the PairInteraction logo using QPainter, ported from data/icon.tex (TikZ source)."""
    s = size / 19.0

    def tx(x: float) -> float:
        return (x + 2.0) * s

    def ty(y: float) -> float:
        return (16.0 - y) * s

    painter.save()
    painter.setRenderHint(QPainter.RenderHint.Antialiasing)
    painter.setClipRect(QRectF(0.0, 0.0, float(size), float(size)))

    # White rounded background
    painter.setPen(Qt.PenStyle.NoPen)
    painter.setBrush(QColor("white"))
    painter.drawRoundedRect(QRectF(0.0, 0.0, float(size), float(size)), 4.0 * s, 4.0 * s)

    # Thick round-capped pen
    pen = QPen()
    pen.setWidthF(0.8 * s)
    pen.setCapStyle(Qt.PenCapStyle.RoundCap)
    pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
    painter.setBrush(Qt.BrushStyle.NoBrush)

    # Coordinate axes
    pen.setColor(QColor("black"))
    painter.setPen(pen)
    axes = QPainterPath()
    axes.moveTo(tx(0), ty(13))
    axes.lineTo(tx(0), ty(0))
    axes.lineTo(tx(15), ty(0))
    painter.drawPath(axes)

    # Red curve
    pen.setColor(QColor("red"))
    painter.setPen(pen)
    red = QPainterPath()
    red.moveTo(tx(1), ty(12))
    red.cubicTo(tx(3), ty(-10), tx(6), ty(2), tx(15), ty(1))
    painter.drawPath(red)

    # Blue curve
    pen.setColor(QColor("blue"))
    painter.setPen(pen)
    blue = QPainterPath()
    blue.moveTo(tx(2), ty(12))
    blue.cubicTo(tx(4), ty(-6), tx(6), ty(5), tx(15), ty(2))
    painter.drawPath(blue)

    # PI text
    font = QFont("serif")
    font.setPixelSize(round(10.0 * s))
    painter.setFont(font)
    pen.setColor(QColor("black"))
    painter.setPen(pen)
    cx, cy = tx(9), ty(9)
    painter.drawText(
        QRectF(cx - 7.0 * s, cy - 6.0 * s, 14.0 * s, 12.0 * s),
        Qt.AlignmentFlag.AlignCenter,
        "PI",
    )

    painter.restore()
