import os, numpy as np
from cdl.encoding import load, parse_aa
HERE = os.path.dirname(__file__)

def test_parse_aa():
    assert parse_aa("D40A:G41S") == {40: "A", 41: "S"}

def test_mutation_list_groups_by_position():
    X, y, names, groups, nom = load(os.path.join(HERE,"fixtures/mini_mutation.csv"), "mutation_list")
    assert nom == len(set(groups))
    assert set(groups) <= {"p40", "p41"}          # one group per mutated position
    assert X.shape[0] == 4 and len(y) == 4

def test_design_matrix_numeric_and_categorical():
    X, y, names, groups, nom = load(os.path.join(HERE,"fixtures/mini_designmatrix.csv"), "design_matrix")
    assert "geneA" in groups                        # numeric col -> its own group
    assert "promoter" in groups                     # categorical -> one group, one-hot minus ref
    assert nom == 2
