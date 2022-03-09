"""
Unit and regression test for the gb_network package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import gb_network


def test_gb_network_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "gb_network" in sys.modules
