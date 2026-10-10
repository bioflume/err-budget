/* Work around the Ubuntu 24.04 Quantum ESPRESSO 6.7 crash ("buffer overflow
   detected" in md5_from_file): the fortified __snprintf_chk aborts on an
   off-by-one terminator write that unfortified builds of QE perform silently.
   Preloading this restores the unfortified behaviour. */
#include <stdarg.h>
#include <stdio.h>
int __snprintf_chk(char *s, size_t maxlen, int flag, size_t slen, const char *fmt, ...) {
    va_list ap; va_start(ap, fmt);
    int r = vsnprintf(s, maxlen, fmt, ap);
    va_end(ap);
    return r;
}
