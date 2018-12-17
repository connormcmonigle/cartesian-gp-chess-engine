# cartesian-gp-chess-engine

This project is an experiment involving evolving cartesian genetic programs to evaluate chess positions. It uses self play to evolve slightly superior chess playing cartesian programs each generation. For each generation, a population is generated through mutations applied to the prior champion. Each mutated cartesian genetic program is evaluated by its performance against the former champion. The highest performing cartesian program is selected for the next round and is deemed the new champion. For each match, moves are selected through a negamax search at depth 2. Despite early signs of promise in which the programs learned simple mates to value queens, the project didn't succeed with cartesian genetic programs proving poor at position evaluation and too slow to be practical.

# Usage

For those interested in training cartesian genetic programs, the chessga.cpp can be compiled with "g++ chessga.cpp -o chessga --std=c++17 -O4". Every 100 generations, a trained cartesian program is saved in a user specified text file. Additionally, in the save.txt file is a somewhat evolved example cartesian program that can be loaded for further evolution. Every 100 generations, the current matches are printed to the console with positions displayed as follows:

eval :: 0.03125
 . . . ♔ . . . ♖
 . ♙ . . . ♙ ♙ .
 . . . ♖ . . . .
 ♙ . ♜ . ♙ ♗ ♟ ♙
 ♕ . . . . . . ♟
 ♝ . . ♙ . . . .
 ♚ . . ♟ . ♛ . .
 . ♞ . . . ♝ ♞ .
 
 Different board representations and cartesian genetic program functions should probably be experimented with. Currently, the cartesian program is fed 25 8x8 board states corresponding to the position. 

# Work In Progess
- Move generation is probably bugged and has not been thoroughly tested. Thanks to Pradyumna Kannan for the magic bitboard implementation.
- Overall code quality is poor due to this project being purely experimental. Many components were hacked on after the initial project was finished.
