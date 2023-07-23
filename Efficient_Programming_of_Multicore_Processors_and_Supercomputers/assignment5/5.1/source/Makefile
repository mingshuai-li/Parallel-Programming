CC = mpicc
CFLAGS = -Wall -Wextra
LDFLAGS = -lm

TARGET = pingpong

SRCS = pingpong.c
OBJS = $(SRCS:.c=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(TARGET) $(LDFLAGS)

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET) job.ou*
