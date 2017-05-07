#ifndef _RANDOM_
#define _RANDOM_ 1

inline
int random_int(const int min, const int max)
{
    int range = max - min;
    return min + int((double(rand())/double(RAND_MAX)) * range);
}

#endif
