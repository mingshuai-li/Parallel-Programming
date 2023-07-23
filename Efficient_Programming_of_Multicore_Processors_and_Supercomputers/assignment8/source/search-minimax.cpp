#include "search.h"
#include "board.h"
#include "eval.h"
#include <stdio.h>
#include <omp.h>
#include <cstring>
#include <chrono>
#include <iostream>

/**
 * To create your own search strategy:
 * - Copy this file into another one.
 * - Change the class name to the name given in the constructor.
 * - Adjust clone() to return an instance of your class.
 * - Adjust the last line of this file to create a global instance of your class.
 * - Adjust the Makefile to include your class in SEARCH_OBJS.
 * - Implement searchBestMove().
 *
 * Advice for implementing searchBestMove():
 * - Call foundBestMove() when finding the best move since the search start.
 * - Call finishedNode() when finishing the evaluation of a tree node.
 * - Use _maxDepth for the strength level (maximal level searched in the tree).
 */
class MinimaxStrategy : public SearchStrategy
{
public:
    // Defines the name of the strategy
    MinimaxStrategy() : SearchStrategy("Minimax") {}

    // Factory method: just return a new instance of this class
    SearchStrategy *clone() { return new MinimaxStrategy(); }

private:
    /**
     * Implementation of the strategy.
     */
    void searchBestMove();

    // Recursive minimax function with alpha-beta pruning
    int minimax(char depth, Board *minimax_board, Evaluator *temp_eval, int alpha, int beta);

    // Recursive minimax function with OpenMP parallelization
    int minimax_omp(char depth, Board minimax_board);

    // Helper function to check if two game boards are the same
    bool isSameBoard(const int* temp_field, const int before_field[121]);

    int adapt_search_depth{5};  // Depth of the search adapted based on game-specific conditions
    int bestEval{-20000};       // Stores the best evaluation value found during the search
    int secondLastField[121];   // Game board field of the second last move
    int previousField[121];     // Game board field of the previous move
    int num_move{1};            // Number of moves played
    Variation _pv;              // Stores the principal variation (best moves sequence)
    Move secondLastMove;        // Second last move played
    Move previousMove;          // Previous move played
    double time_remain{60};     // Remaining time for the search
};

bool MinimaxStrategy::isSameBoard(const int* temp_field, const int before_field[121])
{
    // Check if two game boards are the same
    return std::memcmp(temp_field, before_field, sizeof(int) * 121) == 0;
}

int MinimaxStrategy::minimax(char depth, Board* minimax_board, Evaluator* temp_eval, int alpha, int beta)
{
    bool MaxorMin = !(depth % 2);

    // If the maximum search depth is reached, calculate the evaluation directly
    if (depth >= adapt_search_depth)
    {
        int evalMultiplier = 1 - MaxorMin * 2;
        int eval = evalMultiplier * temp_eval->calcEvaluation(minimax_board);
        return eval;
    }

    int eval;
    MoveList list;
    Move possibleMove;

    // Generate all possible moves for the current board position
    minimax_board->generateMoves(list);

    int current_value = MaxorMin ? -20000 : 20000;

    // Iterate over each move and evaluate the resulting board positions recursively
    int movelength = list.getLength();

    for (int i = movelength - 1; i >= 0; --i)
    {
        list.access(possibleMove, i);
        minimax_board->playMove(possibleMove);
        eval = minimax(depth + 1, minimax_board, temp_eval, alpha, beta);
        minimax_board->takeBack();

        if (MaxorMin)
        {
            current_value = std::max(current_value, eval);
            alpha = std::max(alpha, eval);
            if (current_value >= beta)
                break;
        }
        else
        {
            current_value = std::min(current_value, eval);
            beta = std::min(beta, eval);
            if (current_value <= alpha)
                break;
        }
    }

    return current_value;
}

int MinimaxStrategy::minimax_omp(char depth, Board minimax_board)
{
    int bestEval = -20000;

    MoveList list;

    // Generate all possible moves for the current board position
    minimax_board.generateMoves(list);
    int movelength = list.getLength();

    // Parallelize the move evaluation using OpenMP
    #pragma omp parallel for schedule(dynamic, 1) shared(bestEval) firstprivate(minimax_board)
    for (int i = movelength - 1; i >= 0; --i)
    {
        // Move possibleMove = moves[i];
        Move possibleMove;
        list.access(possibleMove, i);
        int eval;
        Evaluator temp_eval;
        
        minimax_board.playMove(possibleMove);
        eval = minimax(depth + 1, &minimax_board, &temp_eval, -20000, 20000);
        minimax_board.takeBack();
        
        if (eval > bestEval)
        {
            // Update the best evaluation and best move in a critical section
            #pragma omp critical
            {
                bestEval = eval;
                _bestMove = possibleMove;
            }
        }
    }

    return bestEval;
}


void MinimaxStrategy::searchBestMove()
{
    auto starttime = std::chrono::high_resolution_clock::now();

    // Check if the same board position occurred in the previous two moves
    
    if (num_move >= 3 && isSameBoard(_board->fieldArray(), secondLastField)) {
        _bestMove = secondLastMove;
        return;
    } 
    // Choose the adapted search depth based on game-specific conditions
    adapt_search_depth = (time_remain < 4) ? 3 :
                         (time_remain < 10) ? 4 :
                         (time_remain > 10 && bestEval < -400 || (time_remain > 20 && bestEval < -100)) ? 6 :
                         5;
    if (num_move == 1){
        adapt_search_depth = 3;
    }
    // Set the number of OpenMP threads to be used
    omp_set_dynamic(0);
    omp_set_num_threads(48);

    // Perform the minimax search
    bestEval = minimax_omp(0, *_board);

    // Output the search results
    std::cout << "Best Eval = " << bestEval << std::endl;

    auto endtime = std::chrono::high_resolution_clock::now();

    double searchtime = std::chrono::duration<double, std::micro>(endtime - starttime).count(); 
    time_remain -= (searchtime / 1000000.0); 

    std::cout << "Depth = " << adapt_search_depth << ", Remain time = " << time_remain << ", MoveNumber = " << num_move << std::endl;

    if (num_move >= 2) {
        // Update the previous and second-last board fields and moves
        std::memcpy(secondLastField, previousField, 121 * sizeof(int));
        std::memcpy(previousField, _board->fieldArray(), 121 * sizeof(int));
        secondLastMove = previousMove;
        previousMove = _bestMove;
    }

    num_move++;
}


// Register ourselves as a search strategy
MinimaxStrategy minimaxStrategy;
