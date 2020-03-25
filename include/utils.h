#pragma once

#include <cstdlib>

#define ASSERT(expression)                                                     \
    do {                                                                       \
        if (!expression)                                                       \
            std::abort();                                                      \
    } while (false)