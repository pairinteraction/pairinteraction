# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import logging
from functools import partial
from typing import Optional

from PySide6.QtCore import QObject, QTimer
from PySide6.QtGui import QStatusTipEvent
from PySide6.QtWidgets import QApplication


def show_status_tip(parent: QObject, message: str, timeout: int = 0, logger: Optional[logging.Logger] = None) -> None:
    """Show a status tip message using QStatusTipEvent."""
    QApplication.sendEvent(parent, QStatusTipEvent(message))
    if logger:
        logger.info(message)
    if timeout > 0:
        QTimer.singleShot(timeout, partial(reset_status_tip, parent))


def reset_status_tip(parent: QObject) -> None:
    """Hide the status tip by sending an empty message."""
    QApplication.sendEvent(parent, QStatusTipEvent("Ready"))
