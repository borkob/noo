cmake_minimum_required(VERSION 3.5)

project(labs_neighborhood_search LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(labs_neighborhood_search main.cpp)

include(CTest)
set(TARGET_E 1 1 1 1 2 2 7 3 8 12 13 5 10 6 19 15 24 32 25 29 26)
set(TARGET_PSL 1 1 1 1 1 1 2 1 2 2 2 1 2 1 2 2 2 2 2 2 2)
foreach(L RANGE 4 20)
    list(GET TARGET_E ${L} target_e)
    list(GET TARGET_PSL ${L} target_psl)
    add_test(NAME Test${L} COMMAND bash -c "./labs_neighborhood_search 42 10000000 ${L} > out_${L}.txt")
    add_test(NAME Test_E_${L} COMMAND bash -c "grep -q '^E: ${target_e}' out_${L}.txt")
    add_test(NAME Test_PSL_${L} COMMAND bash -c "grep -q '^PSL: ${target_psl}' out_${L}.txt")
endforeach()

set(L 513)
foreach(seed RANGE 1 25)
    add_test(NAME L${L}_${seed} COMMAND bash -c "./labs_neighborhood_search ${seed} 15000000 ${L} > out_${L}_${seed}.txt")
endforeach()
