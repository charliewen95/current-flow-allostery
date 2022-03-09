"""
Unit and regression test for the flowNetwork package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import flowNetwork


def test_flowNetwork_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "flowNetwork" in sys.modules
