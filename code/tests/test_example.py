"""
Unit tests for the example library
"""

from example import example


class TestExample:
    def test_add_numbers(self):
        assert example.add_one(2) == 3
        assert example.add_one(-1) == 0

    def test_multiply_by_two(self):
        assert example.multiply_by_two(2) == 4

    # def test_multiplication(self):
    #     assert 100 == example.multiply(10, 10)
    def test_cbmpy_example_model(self):
        assert example.cbmpy_example_model() == 1
