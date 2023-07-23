/**
 * Send a start position into a game communication channel,
 * and observe positions sent in the channel
 *
 * (C) 2005-2015, Josef Weidendorfer, GPLv2+
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "board.h"
#include "network.h"

/* Global, static vars */
static NetworkLoop l;
static Board myBoard;

/* Start game? */
static bool sendStart = true;

static bool startX = false;
static bool delayedStart = false;

/* Observe game play? */
static bool observeGame = true;

/* only show reachable positions with valid moves */
static bool onlyReachable = false;

/* time limit in seconds (0: no time limit, -1: not set) */
static int secsToPlay = -1;

/* to set verbosity of NetworkLoop implementation */
extern int verbose;

#define DEFAULT_DOMAIN_PORT 23412

/* remote channel */
static char* host = 0;       /* not used on default */
static int rport = DEFAULT_DOMAIN_PORT;

/* local channel */
static int lport = DEFAULT_DOMAIN_PORT;

/* Where to read position to broadcast from? (0: start position) */
static FILE* file = 0;

class MyDomain: public NetworkDomain
{
public:
    MyDomain(int p) : NetworkDomain(p) { sent = 0; }

    void sendBoard(Board*);

protected:
    void received(char* str);
    void newConnection(Connection*);

private:
    Board* sent;
};

void MyDomain::sendBoard(Board* b)
{
    if (b) {
	static char tmp[500];
	sprintf(tmp, "pos %s\n", b->getState());
	printf("%s", tmp+4);
	int state = b->validState();
	printf("%s\n", Board::stateDescription(state));
	broadcast(tmp);
    }
    sent = b;
}

void MyDomain::received(char* str)
{
    if (strncmp(str, "quit", 4)==0) {
	l.exit();
	return;
    }

    if (strncmp(str, "pos ", 4)!=0) return;

    // on receiving remote position, do not broadcast own board any longer
    sent = 0;
    // stop after receiving remote position if we are not observing game play
    if (!observeGame) l.exit();

    if (myBoard.validState() != Board::empty) {

	Board newBoard;
	newBoard.setState(str+4);

	Move m = myBoard.moveToReach(&newBoard, false);
	if (m.type == Move::none) {
	    printf("WARNING: Got a board which is not reachable via a valid move !?\n");
	    if (onlyReachable) return;
	}
	else {
	    if (myBoard.actColor() == Board::color1)
		printf("O draws move '%s'...\n", m.name());
	    else
		printf("X draws move '%s'...\n", m.name());
	}
    }

    myBoard.setState(str+4);
    printf("%s", str+4);
    int state = myBoard.validState();
    printf("%s\n", Board::stateDescription(state));

    switch(state) {
	case Board::timeout1:
	case Board::timeout2:
	case Board::win1:
	case Board::win2:
	    l.exit();
	default:
	    break;
    }
}

void MyDomain::newConnection(Connection* c)
{
    NetworkDomain::newConnection(c);

    if (sent) {
	static char tmp[500];
	int len = sprintf(tmp, "pos %s\n", sent->getState());
	c->sendString(tmp, len);
    }
}

class SendOnEmptyTimer: public NetworkTimer
{
public:
    SendOnEmptyTimer(MyDomain* d, Board* b)
	: NetworkTimer(500) { domain = d; toSend = b; }

protected:
    void timeout(NetworkLoop*);

private:
    MyDomain* domain;
    Board* toSend;
};

void SendOnEmptyTimer::timeout(NetworkLoop*)
{
    if (!myBoard.isValid()) {
	myBoard = *toSend;
	domain->sendBoard(&myBoard);
    }
}


static void printHelp(char* prg, bool printHeader)
{
    if (printHeader)
	printf("Start V 0.3\n"
	       "Broadcast a game position, observe the game and quit in winning state.\n\n");
    
    printf("Usage: %s [options] [O|X|<file>|-]\n\n"
	   "  O                Use regular start position, O to play (default)\n"
	   "  X                Use regular start position, X to play\n"
	   "  <file>           File containing start position\n"
	   "  -                Position is read from standard input\n\n",
	   prg);
    printf(" Options:\n"
	   "  -h / --help      Print this help text\n"
	   "  -v / -vv         Be verbose / more verbose\n"
	   "  -o               Only observe, no start\n"
	   "  -d               Only start if not receiving other board after 500ms\n"
	   "  -n               Do not observe, but stop when other board received\n"
	   "  -t <timeToPlay>  Start in tournament modus (limited time)\n"
	   "  -r               Only accept positions reachable (default: all)\n"
	   "  -p [host:][port] Connection to broadcast channel\n"
	   "                   (default: %d)\n\n", DEFAULT_DOMAIN_PORT);
    exit(1);
}

static void parseArgs(int argc, char* argv[])
{
    verbose = 0;
    int arg=0;
    while(arg+1<argc) {
	arg++;
	if (argv[arg][0] == '-') {
	    if (strcmp(argv[arg],"-h")==0 ||
		strcmp(argv[arg],"--help")==0) printHelp(argv[0], true);
	    if (strcmp(argv[arg],"-v")==0) {
		verbose++;
		continue;
	    }
	    if (strcmp(argv[arg],"-vv")==0) {
		verbose += 2;
		continue;
	    }
	    if (strcmp(argv[arg],"-o")==0) {
		sendStart = false;
		continue;
	    }
	    if (strcmp(argv[arg],"-n")==0) {
		observeGame = false;
		continue;
	    }
	    if (strcmp(argv[arg],"-d")==0) {
		delayedStart = true;
		continue;
	    }
	    if (strcmp(argv[arg],"-r")==0) {
		onlyReachable = true;
		continue;
	    }
	    if ((strcmp(argv[arg],"-t")==0) && (arg+1<argc)) {
		arg++;
		secsToPlay = atoi(argv[arg]);
		if (secsToPlay == 0) {
		    printf("%s: WARNING - Ignoring tournament; %d secs to play\n",
			   argv[0], secsToPlay);
		}
		continue;
	    }
	    if ((strcmp(argv[arg],"-p")==0) && (arg+1<argc)) {
		arg++;
		if (argv[arg][0]>'0' && argv[arg][0]<='9') {
		    lport = atoi(argv[arg]);
		    continue;
		}
		char* c = strrchr(argv[arg],':');
		int p = 0;
		if (c != 0) {
		    *c = 0;
		    p = atoi(c+1);
		}
		host = argv[arg];
		if (p) rport = p;
		continue;
	    }
	    if (strcmp(argv[arg],"-")==0) {
		file = stdin;
		continue;
	    }
	    printf("%s: ERROR - Unknown option %s\n", argv[0], argv[arg]);
	    printHelp(argv[0], false);
	}

	if (strcmp(argv[arg],"X")==0) {
	    startX = true;
	    continue;
	}
	if (strcmp(argv[arg],"O")==0) {
	    startX = false;
	    continue;
	}

	if (!file) {
	    file = fopen(argv[arg], "r");
	    if (!file) {
		printf("%s: ERROR - Can not open '%s' for reading start position\n",
		       argv[0], argv[arg]);
		printHelp(argv[0], false);
	    }
	}
	break;
    }
    if ((sendStart == false) && (observeGame == false)) {
	printf("Nothing to do, terminating (no board to send and not observing)\n");
	exit(1);
    }
}

int main(int argc, char* argv[])
{
    parseArgs(argc, argv);
    myBoard.setVerbose(verbose);

    MyDomain d(lport);
    l.install(&d);

    if (host) d.addConnection(host, rport);

    if (sendStart) {
	static Board boardToSend;
	if (file) {
	    char tmp[500];
	    int len = 0, c;
	    while( len<499 && (c=fgetc(file)) != EOF)
		tmp[len++] = (char) c;
	    tmp[len++]=0;

	    if (!boardToSend.setState(tmp)) {
		printf("%s: WARNING - Can not parse given position; using start position\n", argv[0]);
		boardToSend.begin(Board::color1);
	    }
	}
	else
	    boardToSend.begin(startX ? Board::color2 : Board::color1);

	if (secsToPlay >= 0) {
	    boardToSend.setMSecsToPlay(Board::color1, 1000 * secsToPlay);
	    boardToSend.setMSecsToPlay(Board::color2, 1000 * secsToPlay);
	}
	if (delayedStart)
	    l.install(new SendOnEmptyTimer(&d, &boardToSend));
	else {
	    myBoard = boardToSend;
	    d.sendBoard(&myBoard);
	}
    }
    return l.run();
}
