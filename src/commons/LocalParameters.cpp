#include "LocalParameters.h"

#ifdef __APPLE__
#include <mach-o/dyld.h>
#include <stdlib.h>
#else
#include <unistd.h>
#endif

std::string LocalParameters::getAbsExePath() {
#ifdef __APPLE__
    uint32_t bufSize = 0;
    // determine required buffer size.
    _NSGetExecutablePath(NULL, &bufSize);
    char* pathBuffer = new char[bufSize];
    if (_NSGetExecutablePath(pathBuffer, &bufSize) != 0) {
        delete[] pathBuffer;
        return std::string();
    }
    char* resolved = realpath(pathBuffer, NULL);
    std::string result = resolved ? resolved : pathBuffer;
    free(resolved);
    delete[] pathBuffer;
    return result;
#else
    const char* linkPath = "/proc/self/exe";
    size_t bufSize = PATH_MAX;
    while (true) {
        std::vector<char> buf(bufSize);
        ssize_t len = readlink(linkPath, buf.data(), bufSize);
        if (len < 0) {
            return std::string();
        }
        if (static_cast<size_t>(len) < bufSize) {
            buf[len] = '\0';  // null-terminate
            char* resolved = realpath(buf.data(), NULL);
            std::string result = resolved ? resolved : std::string(buf.data());
            free(resolved);
            return result;
        }
        bufSize *= 2;
    }
#endif
}
