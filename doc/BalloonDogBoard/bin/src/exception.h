/*
Copyright 2012-2018 Tamas Bolner

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#pragma once

#include <cinttypes>

#define FE_MESSAGE_BUFFER_SIZE 1024

/**
 * Exception class with printf like text formatting
 * capabilities, using variable argument list.
 */
class FException {
private:
    char error_msg[FE_MESSAGE_BUFFER_SIZE];
    int64_t error_code;

protected:

public:
    FException(const char* error_msg, ...);
    FException(int64_t error_code, const char* error_msg, ...);

    void Print() const;
    int64_t getErrorCode() const;
    const char *getMessage() const;
};
