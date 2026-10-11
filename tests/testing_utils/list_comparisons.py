from typing import List


def equivalent_lists(l_a: list, l_b: list) -> bool:
    if len(l_a) != len(l_b):
        return False
    for i in range(len(l_a)):
        if l_a[i] != l_b[i]:
            return False
    return True


def list_in_order(l: list, order: List[int]) -> list:
    ret = []
    for o in order:
        ret.append(l[o])
    return ret

