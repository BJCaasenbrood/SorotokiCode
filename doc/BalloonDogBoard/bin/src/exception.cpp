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

#include <cinttypes>
#include <iostream>
#include <cstdarg>
#include "exception.h"

/**
 * Constructor
 * Use the same parameters as for the printf() function.
 *
 * @param error_msg 
 * @param ... Variable argument list.
 * @returns FException
 */
FException::FException(const char *error_msg, ...) {
    va_list arg_list;

    va_start(arg_list, error_msg);
    vsnprintf(this->error_msg, (size_t)FE_MESSAGE_BUFFER_SIZE, error_msg, arg_list);
    va_end(arg_list);

    this->error_code = 0;
}

/**
 * Constructor
 * After the error_code use the same parameters as for the printf() function.
 * 
 * @param error_code
 * @param error_msg 
 * @param ... Variable argument list.
 * @returns FException
 */
FException::FException(int64_t error_code, const char *error_msg, ...) {
    va_list arg_list;

    va_start(arg_list, error_msg);
    vsnprintf(this->error_msg, (size_t)FE_MESSAGE_BUFFER_SIZE, error_msg, arg_list);
    va_end(arg_list);

    this->error_code = error_code;
}

/**
 * Prints out the error message text to the console.
 */
void FException::Print() const {
    std::cout << this->error_msg;
}

/**
 * Returns the error code if it's set. Otherwise returns 0.
 * 
 * @returns int64_t
 */
int64_t FException::getErrorCode() const {
    return this->error_code;
}

/**
 * Returns the error message text of the exception.
 * 
 * @returns char*
 */
const char* FException::getMessage() const {
    return this->error_msg;
}
