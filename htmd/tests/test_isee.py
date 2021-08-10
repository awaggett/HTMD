"""
Unit and regression test for the htmd package.
"""

# Import package, test suite, and other packages as needed
import htmd_
import pytest
import sys

def test_isee_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "htmd" in sys.modules
