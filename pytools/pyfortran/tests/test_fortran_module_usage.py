from fortran_module_usage.check_fortran_module_usage import (
    join_continued_lines,
    find_unused_subroutines,
    find_unused_routines_from_other_modules,
)

def test_join_continued_lines():
    lines = [
        "test1 & !hello",
        "  (arg1, arg2)",
        "test2",
    ]
    expected = [
        "test1    (arg1, arg2)",
        "test2",
    ]
    assert join_continued_lines(lines) == expected

def test_find_unused_subroutines():
    lines = [
        "subroutine test1",
        "end subroutine test1",
        "subroutine test2",
        "call test1",
        "end subroutine test2",
        "subroutine test3",
        "! call test2",
        "module procedure test3",
        "subroutine test4",
        "public :: test4",
    ]
    unused = find_unused_subroutines(lines)
    assert unused == ["test2"]

def test_find_unused_routines_from_other_modules():
    lines = [
        "use some_module, only: routine1, routine2, routine3",
        "call routine1",
        "! call routine2",
        "call routine3",
        "subroutine routine4",
    ]
    unused = find_unused_routines_from_other_modules(lines)
    assert unused == ["routine2"]

def main():
    test_join_continued_lines()
    test_find_unused_subroutines()
    test_find_unused_routines_from_other_modules()
    print("All tests passed.")

if __name__ == "__main__":
    main()