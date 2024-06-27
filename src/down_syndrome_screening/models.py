# src/down_syndrome_screening/models.py

import numpy as np
from pydantic import BaseModel, Field
from dataclasses import dataclass
from typing import List, Tuple


@dataclass
class Marker:
    name: str
    median_mom_down: float
    median_mom_control: float
    log_sd_down: float
    log_sd_control: float


class Population(BaseModel):
    maternal_age_mean: float = Field(27.0, description="Mean maternal age in years")
    maternal_age_sd: float = Field(
        5.5, description="Standard deviation of maternal age in years"
    )
    markers: List[Marker] = Field(
        default_factory=list, description="List of screening markers"
    )
    down_syndrome_prevalence: float = Field(
        1 / 700, description="Prevalence of Down syndrome in the population"
    )
    correlation_matrix_down: np.ndarray = Field(
        default_factory=lambda: np.eye(3),
        description="Correlation matrix for Down syndrome cases",
    )
    correlation_matrix_control: np.ndarray = Field(
        default_factory=lambda: np.eye(3),
        description="Correlation matrix for control cases",
    )

    class Config:
        arbitrary_types_allowed = True

    def generate_sample(self, size: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        is_down = np.random.random(size) < self.down_syndrome_prevalence
        maternal_ages = np.random.normal(
            self.maternal_age_mean, self.maternal_age_sd, size
        )

        mean_down = np.log([m.median_mom_down for m in self.markers])
        mean_control = np.log([m.median_mom_control for m in self.markers])
        cov_down = np.diag([m.log_sd_down**2 for m in self.markers])
        cov_control = np.diag([m.log_sd_control**2 for m in self.markers])

        cov_down = cov_down @ self.correlation_matrix_down @ cov_down
        cov_control = cov_control @ self.correlation_matrix_control @ cov_control

        log_marker_values_down = np.random.multivariate_normal(
            mean_down, cov_down, np.sum(is_down)
        )
        log_marker_values_control = np.random.multivariate_normal(
            mean_control, cov_control, np.sum(~is_down)
        )

        log_marker_values = np.zeros((size, len(self.markers)))
        log_marker_values[is_down] = log_marker_values_down
        log_marker_values[~is_down] = log_marker_values_control

        marker_values = np.exp(log_marker_values)

        return maternal_ages, marker_values, is_down


def create_default_population() -> Population:
    markers = [
        Marker("Free β-hCG", 1.70, 1.01, 0.28, 0.27),
        Marker("PAPP-A", 0.49, 1.00, 0.31, 0.25),
        Marker("NT", 1.74, 1.01, 0.23, 0.13),
    ]

    correlation_matrix_down = np.array(
        [[1.000, 0.191, 0.000], [0.191, 1.000, 0.000], [0.000, 0.000, 1.000]]
    )

    correlation_matrix_control = np.array(
        [[1.000, 0.186, 0.000], [0.186, 1.000, 0.000], [0.000, 0.000, 1.000]]
    )

    return Population(
        markers=markers,
        correlation_matrix_down=correlation_matrix_down,
        correlation_matrix_control=correlation_matrix_control,
    )

