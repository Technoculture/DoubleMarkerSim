# src/down_syndrome_screening/models.py

from dataclasses import dataclass
from typing import List, Tuple
from pydantic import BaseModel, Field
import numpy as np
from scipy.stats import multivariate_normal


@dataclass
class Marker:
    name: str
    median_mom_down: float
    median_mom_control: float
    mean_log_mom_down: float
    mean_log_mom_control: float
    sd_down: float
    sd_control: float


class Population(BaseModel):
    maternal_age_mean: float = Field(27.0, description="Mean maternal age in years")
    maternal_age_sd: float = Field(
        5.5, description="Standard deviation of maternal age in years"
    )
    markers: List[Marker] = Field(
        default_factory=list, description="List of screening markers"
    )
    down_syndrome_prevalence: float = Field(
        float(1 / 700), description="Prevalence of Down syndrome in the population"
    )
    correlation_matrix_down: np.ndarray = Field(
        default_factory=lambda: np.eye(7),
        description="Correlation matrix for Down syndrome cases",
    )
    correlation_matrix_control: np.ndarray = Field(
        default_factory=lambda: np.eye(7),
        description="Correlation matrix for control cases",
    )

    class Config:
        arbitrary_types_allowed = True

    def generate_sample(self, size: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Generate a sample of the population using multivariate log-normal distribution.

        Args:
            size (int): Number of individuals to generate.

        Returns:
            Tuple[np.ndarray, np.ndarray, np.ndarray]: Arrays of maternal ages, marker values, and Down syndrome status.
        """
        # Generate Down syndrome status
        is_down = np.random.random(size) < self.down_syndrome_prevalence

        # Generate maternal ages
        maternal_ages = np.random.normal(
            self.maternal_age_mean, self.maternal_age_sd, size
        )

        # Prepare mean vectors and covariance matrices
        mean_down = np.array([m.mean_log_mom_down for m in self.markers])
        mean_control = np.array([m.mean_log_mom_control for m in self.markers])
        cov_down = np.diag([m.sd_down**2 for m in self.markers])
        cov_control = np.diag([m.sd_control**2 for m in self.markers])

        # Apply correlation matrices
        cov_down = cov_down @ self.correlation_matrix_down @ cov_down
        cov_control = cov_control @ self.correlation_matrix_control @ cov_control

        # Generate marker values using multivariate normal distribution
        marker_values_down = multivariate_normal.rvs(
            mean=mean_down, cov=cov_down, size=np.sum(is_down)
        )
        marker_values_control = multivariate_normal.rvs(
            mean=mean_control, cov=cov_control, size=np.sum(~is_down)
        )

        # Combine Down syndrome and control marker values
        marker_values = np.zeros((size, len(self.markers)))
        marker_values[is_down] = marker_values_down
        marker_values[~is_down] = marker_values_control

        return maternal_ages, marker_values, is_down


def create_default_population() -> Population:
    markers = [
        Marker("PAPP-A", 0.49, 1.00, -0.32, 0.00, 0.31, 0.25),
        Marker("Free β-hCG", 1.70, 1.01, 0.24, 0.01, 0.28, 0.27),
        Marker("ADAM12", 0.89, 1.00, -0.07, -0.01, 0.19, 0.16),
        Marker("Total hCG", 1.28, 1.00, 0.11, -0.02, 0.20, 0.19),
        Marker("PP13", 0.91, 1.00, -0.03, -0.01, 0.17, 0.19),
        Marker("PlGF", 0.80, 1.00, -0.10, -0.02, 0.14, 0.14),
        Marker("NT", 1.74, 1.01, 0.24, 0.00, 0.23, 0.13),
    ]

    # Correlation matrices from the paper (Table 3)
    correlation_matrix_down = np.array(
        [
            [1.000, 0.191, 0.460, 0.182, 0.408, 0.152, 0.000],
            [0.191, 1.000, 0.297, 0.715, 0.389, -0.124, 0.000],
            [0.460, 0.297, 1.000, 0.598, 0.528, 0.154, 0.000],
            [0.182, 0.715, 0.598, 1.000, 0.627, -0.022, 0.000],
            [0.408, 0.389, 0.528, 0.627, 1.000, 0.018, 0.000],
            [0.152, -0.124, 0.154, -0.022, 0.018, 1.000, 0.000],
            [0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 1.000],
        ]
    )

    correlation_matrix_control = np.array(
        [
            [1.000, 0.186, 0.413, 0.221, 0.324, 0.256, 0.000],
            [0.186, 1.000, 0.152, 0.677, 0.287, 0.102, 0.000],
            [0.413, 0.152, 1.000, 0.434, 0.432, 0.324, 0.000],
            [0.221, 0.677, 0.434, 1.000, 0.531, 0.189, 0.000],
            [0.324, 0.287, 0.432, 0.531, 1.000, 0.216, 0.000],
            [0.256, 0.102, 0.324, 0.189, 0.216, 1.000, 0.000],
            [0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 1.000],
        ]
    )

    return Population(
        markers=markers,
        correlation_matrix_down=correlation_matrix_down,
        correlation_matrix_control=correlation_matrix_control,
    )
