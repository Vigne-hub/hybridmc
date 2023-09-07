import sys


def make_rc_tuple(nonlocal_bonds, rc):
    """
    Function to produce nonlocal bonds list with each nonlocal bond list element containing the rc as the third
    element

    args:
    nonlocal_bonds -- List of Lists containing the indices of beads bonded to each other, for example [ [1, 9], [4, 15]].
    These elements may potentially have the bond length (rc) as the third element as well or not,
    for example [[1, 9, 2.3], [2, 6]].

    rc -- Default rc to be appended to each nonlocal_bonds element in case they do not have this information

    returns:
    The same list but with each list element having the bond length (rc) as their third element.

    """

    # loop through each bonded bead pair
    for el in nonlocal_bonds:

        # check if rc for the bond is already given
        if len(el) != 3:
            # if not then append rc
            el.append(rc)

    return nonlocal_bonds


def sort_triplet(bond_list):
    """
swap bond i and j in (i, j, rc) triplets within bond list if i > j

params:

bond_list: List[List]  - bond list

"""

    for el in bond_list:
        if el[0] > el[1]:
            el[0], el[1] = el[1], el[0]

    return bond_list


def stairs_check(bond_pair, stair_rcs, nonlocal_bonds):
    """
    function to check if staircase distances are appropriate

    """
    rc_check = 0
    for el in nonlocal_bonds:
        if el[:-1] == bond_pair:
            rc_check = el[-1]
    problem_rc = []
    for el in stair_rcs:
        if el < rc_check:
            problem_rc.append(el)
    if problem_rc:
        sys.exit((
            f"Value error: rc values given ({problem_rc}) are too small, no rc can be less than smallest allowed ({rc_check})"))


def bits_to_bonds(bits, nonlocal_bonds):
    bits = [i for i, bit in enumerate(bits) if bit]
    bonds = [nonlocal_bonds[i] for i in bits]
    return bonds


def format_bits(bits):
    return ''.join(map(lambda x: '1' if x else '0', bits))
