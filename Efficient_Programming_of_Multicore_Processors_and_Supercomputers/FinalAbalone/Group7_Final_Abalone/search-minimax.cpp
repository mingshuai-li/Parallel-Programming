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
 * - Change the class name to the desired name in the constructor.
 * - Adjust the `clone()` method to return an instance of your class.
 * - Adjust the last line of this file to create a global instance of your class.
 * - Adjust the Makefile to include your class in `SEARCH_OBJS`.
 * - Implement the `searchBestMove()` method.
 *
 * Advice for implementing `searchBestMove()`:
 * - Call `foundBestMove()` when finding the best move since the start of the search.
 * - Call `finishedNode()` when finishing the evaluation of a tree node.
 * - Use `_maxDepth` for the strength level (the maximum level searched in the tree).
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

    // Minimax algorithm with OpenMP parallelization
    int minimax_omp(int depth, Board minimax_board);

    // Minimax algorithm implementation
    int minimax(int depth, Board *minimax_board, Evaluator *temp_eval, int alpha, int beta);

    // Function to check if two boards are the same
    bool isSameBoard(const int *temp_field, const int *before_field);

    int bestEval;                  // Stores the best evaluation value
    int adapt_search_depth{5};     // Adapted search depth based on remaining time and evaluation
    int num_move{1};               // Number of moves made so far
    int secondLastField[121];      // Field configuration before the second last move
    int previousField[121];        // Field configuration before the previous move
    Move secondLastMove;           // Second last move
    Move previousMove;             // Previous move
    double time_remain{60};        // Remaining time in seconds
};

// Function to check if two boards are the same
bool MinimaxStrategy::isSameBoard(const int *temp_field, const int *before_field)
{
    return std::memcmp(temp_field, before_field, sizeof(int) * 121) == 0;
}

// Minimax algorithm implementation
int MinimaxStrategy::minimax(int depth, Board *minimax_board, Evaluator *temp_eval, int alpha, int beta)
{
    bool MaxorMin = !(depth % 2);

    if (depth >= adapt_search_depth)
    {
        int evalMultiplier = 1 - MaxorMin * 2;
        return evalMultiplier * temp_eval->calcEvaluation(minimax_board);
    }

    MoveList list;
    Move possibleMove;
    minimax_board->generateMoves(list);

    int currentValue = MaxorMin ? -20000 : 20000;

    while (list.getNext(possibleMove))
    {
        minimax_board->playMove(possibleMove);
        int eval = minimax(depth + 1, minimax_board, temp_eval, alpha, beta);
        minimax_board->takeBack();

        if (MaxorMin)
        {
            currentValue = std::max(currentValue, eval);
            alpha = std::max(alpha, eval);
            if (currentValue >= beta)
                break;
        }
        else
        {
            currentValue = std::min(currentValue, eval);
            beta = std::min(beta, eval);
            if (currentValue <= alpha)
                break;
        }
    }

    return currentValue;
}

// Minimax algorithm with OpenMP parallelization
int MinimaxStrategy::minimax_omp(int depth, Board minimax_board)
{
    int bestEval = -20000;

    MoveList list;
    
    // Generate all possible moves for the current board position
    minimax_board.generateMoves(list);
    int movelength = list.getLength();

    // Parallelize the search using OpenMP
    #pragma omp parallel for schedule(dynamic, 1) shared(bestEval) firstprivate(minimax_board)
    for (int i = movelength - 1; i >= 0; --i)
    {
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

// Main function for searching the best move
void MinimaxStrategy::searchBestMove()
{
    auto starttime = std::chrono::high_resolution_clock::now();

    if (num_move >= 3 && isSameBoard(_board->fieldArray(), secondLastField))
    {
        // Use the second last move if the board configuration is the same
        _bestMove = secondLastMove;
    }
    else
    {
        // Determine the adapted search depth based on various factors
        adapt_search_depth = (time_remain < 5) ? 3 :
                              (time_remain < 10) ? 4 :
                              (bestEval < -600) ? 6 :
                              (time_remain > 20 && bestEval < -400) ? 6 :
                              5;

        omp_set_dynamic(0);
        omp_set_num_threads(48);
        bestEval = minimax_omp(0, *_board);

        //std::cout << "Final best Eval = " << bestEval << std::endl;
    }

    auto endtime = std::chrono::high_resolution_clock::now();

    double searchtime = std::chrono::duration<double, std::micro>(endtime - starttime).count();

    time_remain -= (searchtime / 1000000.0);

    //std::cout << "AdaptDepth = " << static_cast<int>(adapt_search_depth)
    //      << ", Remaining time = " << time_remain
    //      << "s, OwnMoveNumber = " << num_move << std::endl;

    if (num_move >= 2)
    {
        // Update the board configurations and move history
        std::memcpy(secondLastField, previousField, 121 * sizeof(int));
        std::memcpy(previousField, _board->fieldArray(), 121 * sizeof(int));
        secondLastMove = previousMove;
        previousMove = _bestMove;
    }

    num_move++;
}

// Register ourselves as a search strategy
MinimaxStrategy minimaxStrategy;

