# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtGui import QFont, QPixmap
from PySide6.QtWidgets import QFrame, QLabel, QScrollArea, QVBoxLayout, QWidget

from pairinteraction_gui.page.base_page import BasePage

images_dir = Path(__file__).parent.parent / "images"


class AboutPage(BasePage):
    """Page for displaying information about pairinteraction and the pairinteraction gui."""

    title = "About"
    tooltip = "Learn more about pairinteraction"
    # icon_path = Path(__file__).parent.parent / "icons" / "about.svg"

    def setupWidget(self) -> None:
        """Set up the about page with application information."""
        # Create a scroll area to handle content that might overflow
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)

        # Create a container widget for the scroll area
        container = QWidget()
        layout = QVBoxLayout(container)
        layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        layout.setSpacing(30)
        layout.setContentsMargins(40, 40, 40, 40)

        # Title with modern styling
        title = QLabel("Pairinteraction")
        title_font = QFont()
        title_font.setPointSize(32)
        title_font.setBold(True)
        title.setFont(title_font)
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(title)

        # Description with better formatting
        desc = QLabel(
            "<p style='font-size: 14px; line-height: 1.5;'>"
            "The pairinteraction software calculates properties of Rydberg atoms. "
            "It consists of a Python library and this graphical user interface for obtaining "
            "single-atom properties and calculating pair potentials, making use of a "
            "high-performance C++ backend."
            "</p>"
        )
        desc.setWordWrap(True)
        desc.setTextFormat(Qt.TextFormat.RichText)
        layout.addWidget(desc)

        # Features with improved styling
        features = QLabel(
            "<div style='font-size: 14px; line-height: 1.6;'>"
            "<b>Key Features:</b><br>"
            "• Optimized construction and diagonalization of Hamiltonians<br>"
            "• Support for single-channel quantum defect theory (SQDT)<br>"
            "• Support for multi-channel quantum defect theory (MQDT)<br>"
            "• Electric and magnetic fields in arbitrary directions<br>"
            "• Support for diamagnetism"
            "</div>"
        )
        features.setWordWrap(True)
        features.setTextFormat(Qt.TextFormat.RichText)
        layout.addWidget(features)

        # Benchmark image
        image_label = QLabel()
        pixmap = QPixmap(images_dir / "0845d67063_1.4.2-cp313-win_amd-ryzen-7-5700g-with-radeon-graphics_reps4.png")
        scaled_pixmap = pixmap.scaled(
            800, 400, Qt.AspectRatioMode.KeepAspectRatio, Qt.TransformationMode.SmoothTransformation
        )
        image_label.setPixmap(scaled_pixmap)
        image_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(image_label)

        # Image caption
        caption = QLabel(
            "<p style='font-size: 12px; color: #666; text-align: center;'>"
            "Figure: Benchmarking the construction and diagonalization of a Hamiltonian of a pair of Rb 60S atoms "
            "for 100 different internuclear distances on an AMD Ryzen 7 5700G CPU using Windows 11. "
            "The Hilbert space comprises pair states that differ at most by 4 in n, l and 25GHz in energy. "
            "When supported, symmetries where used to reduce the Hilbert space size. "
            "See the <a href='https://github.com/pairinteraction/pairinteraction/tree/master/tools/benchmarking' "
            "style='color: #007bff;'>benchmarking tool</a>."
            "</p>"
        )
        caption.setWordWrap(True)
        caption.setTextFormat(Qt.TextFormat.RichText)
        caption.setAlignment(Qt.AlignmentFlag.AlignCenter)
        caption.setOpenExternalLinks(True)
        layout.addWidget(caption)

        # Citation with better styling
        citation = QLabel(
            "<div style='font-size: 14px; line-height: 1.5; margin: 20px 0;'>"
            "<b>Please Cite:</b><br>"
            "Sebastian Weber, Christoph Tresp, Henri Menke, Alban Urvoy, Ofer Firstenberg, "
            "Hans Peter Büchler, Sebastian Hofferberth, "
            "<i>Tutorial: Calculation of Rydberg interaction potentials</i>, "
            "J. Phys. B: At. Mol. Opt. Phys. 50, 133001 (2017)"
            "</div>"
        )
        citation.setWordWrap(True)
        citation.setTextFormat(Qt.TextFormat.RichText)
        citation.setOpenExternalLinks(True)
        layout.addWidget(citation)

        # Links with improved styling
        links = QLabel(
            "<div style='font-size: 14px; line-height: 1.5;'>"
            "<b>Links:</b><br>"
            "• <a href='https://github.com/pairinteraction/pairinteraction' "
            "style='color: #007bff;'>GitHub Repository</a><br>"
            "• <a href='https://arxiv.org/abs/1612.08053' style='color: #007bff;'>arXiv:1612.08053</a><br>"
            "• <a href='https://www.gnu.org/licenses/lgpl-3.0.html' style='color: #007bff;'>License: LGPL v3</a>"
            "</div>"
        )
        links.setTextFormat(Qt.TextFormat.RichText)
        links.setOpenExternalLinks(True)
        layout.addWidget(links)

        # Set the container as the scroll area widget
        scroll.setWidget(container)

        # Add the scroll area to the main layout
        self.layout().addWidget(scroll)
