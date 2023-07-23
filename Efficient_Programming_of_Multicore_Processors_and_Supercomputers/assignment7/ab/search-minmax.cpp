#include "search.h"
#include "board.h"
#include "eval.h"
#include <stdio.h>
#include <omp.h>
#include <iostream>
#include <chrono>

/**
 * To create your own search strategy:
 * - copy this file into another one,
 * - change the class name one the name given in constructor,
 * - adjust clone() to return an instance of your class
 * - adjust last line of this file to create a global instance
 *   of your class
 * - adjust the Makefile to include your class in SEARCH_OBJS
 * - implement searchBestMove()
 *
 * Advises for implementation of searchBestMove():
 * - call foundBestMove() when finding a best move since search start
 * - call finishedNode() when finishing evaluation of a tree node
 * - Use _maxDepth for strength level (maximal level searched in tree)
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
    /* recursive minimax search top layer*/
    int minimax_omp(int depth, Board minimax_board, int& num_eval);
    /* recursive minimax search */
    int minimax(int depth, Board * minimax_board, Evaluator * ev, int& num_eval);

    int _lastBestEval;

    Variation _pv;

    int maxdepth;
    short _ownMoveNumber{1};
};

int MinimaxStrategy::minimax(int depth, Board* minimax_board, Evaluator* ev, int& num_eval)
{
    bool maximizeTurn = !(depth % 2); 

    if (depth >= maxdepth) 
    {
        int eval = (1 - maximizeTurn * 2) * ev->calcEvaluation(minimax_board);
        num_eval++;
        return eval;
    }

    int eval;
    MoveList list;
    Move m;

    minimax_board->generateMoves(list);

    if(maximizeTurn){
        int bestValue = -20000; 
        while(list.getNext(m))
        {
            minimax_board->playMove(m);
            eval = minimax(depth + 1, minimax_board, ev, num_eval);
            minimax_board->takeBack();
            if(eval > bestValue){
                bestValue = eval;
            }
        }
        return bestValue;
    }
    else{
        int worstValue = 20000; 
        while(list.getNext(m))
        {
            minimax_board->playMove(m);
            eval = minimax(depth + 1, minimax_board, ev, num_eval);
            minimax_board->takeBack();

            if(eval < worstValue){
                worstValue = eval;
            }
        }
        return worstValue;
    }
}

int MinimaxStrategy::minimax_omp(int depth, Board minimax_board, int& num_eval)
{
    bool maximizeTurn = true; 
    int bestEval = -20000;

    MoveList list;
    Move moves[150];

    minimax_board.generateMoves(list);
    int num_moves = list.getLength();

    for(int i = 0; i < num_moves; ++i){
        list.getNext(moves[i]);
    }

    #pragma omp parallel for schedule(dynamic, 1) reduction(+: num_eval) shared(bestEval) firstprivate(minimax_board)
    for(int i = 0; i < num_moves; ++i)
    {
        Move m = moves[i];
        int eval;
        Evaluator ev;

        minimax_board.playMove(m);
        eval = minimax(depth + 1, &minimax_board, &ev, num_eval);
        minimax_board.takeBack();

        if (eval > bestEval)
        {
            #pragma omp critical
            {
            bestEval = eval;
            _bestMove = m;
            }
        }
    }

    return bestEval;
}

void MinimaxStrategy::searchBestMove()
{
    int num_eval = 0;
    Move m;

    auto t1 = std::chrono::high_resolution_clock::now();

    maxdepth = _maxDepth;

    omp_set_dynamic(0);
    omp_set_num_threads(2);
    
    _lastBestEval = minimax_omp(0, *_board, num_eval); 

    std::cout << " Best Eval = " << _lastBestEval << std::endl;
    std::cout << "Number of Evaluations = " << num_eval << std::endl;

    auto t2 = std::chrono::high_resolution_clock::now();

    double searchtime = std::chrono::duration<double, std::micro>(t2 - t1).count();
    
    double evalPerSecond = num_eval / searchtime;

    std::cout << "Evaluations per second = " << evalPerSecond << " * 10^6" << std::endl;

    _ownMoveNumber++;
}


// register ourselve as a search strategy
MinimaxStrategy minimaxStrategy;