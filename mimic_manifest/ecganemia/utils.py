"""Ortak yardımcılar: seed, logging, zamanlama. (identity/ecgid/utils.py ile aynı.)"""
from __future__ import annotations

import logging
import random
import time
from contextlib import contextmanager

import numpy as np


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    try:
        import torch

        torch.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    except ImportError:
        pass


def get_logger(name: str) -> logging.Logger:
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(
            logging.Formatter("%(asctime)s | %(levelname)s | %(name)s | %(message)s",
                              datefmt="%H:%M:%S")
        )
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
    return logger


@contextmanager
def timer(logger: logging.Logger, label: str):
    t0 = time.perf_counter()
    logger.info("BAŞLADI: %s", label)
    yield
    logger.info("BİTTİ: %s (%.1f s)", label, time.perf_counter() - t0)
