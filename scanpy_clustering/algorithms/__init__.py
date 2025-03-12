"""
Algorithm registry for modular clustering implementations
"""
from typing import Dict, Type

from .base import BaseAlgorithm

# Algorithm registry to be populated
_ALGORITHMS: Dict[str, Type[BaseAlgorithm]] = {}

def register_algorithm(name: str, algorithm_class: Type[BaseAlgorithm]) -> None:
    """
    Register a new algorithm implementation.
    
    Parameters
    ----------
    name : str
        Name of the algorithm.
    algorithm_class : Type[BaseAlgorithm]
        Algorithm class.
    """
    _ALGORITHMS[name] = algorithm_class

def get_algorithm(name: str) -> BaseAlgorithm:
    """
    Get algorithm implementation by name.
    
    Parameters
    ----------
    name : str
        Name of the algorithm.
        
    Returns
    -------
    BaseAlgorithm
        Algorithm implementation.
        
    Raises
    ------
    ValueError
        If algorithm is not registered.
    """
    if name not in _ALGORITHMS:
        if len(_ALGORITHMS) == 0:
            raise ValueError(
                f"No algorithms registered yet. Please implement and register an algorithm."
            )
        raise ValueError(
            f"Unknown algorithm: {name}. "
            f"Available algorithms: {list(_ALGORITHMS.keys())}"
        )
    return _ALGORITHMS[name]() 