/**
 * @file socket_api.h
 * @version 0.1
 * @date 2024-08-10
 * @copyright Copyright (c) 2024
 * 
 */


#pragma once
#include <iostream>
#include <cstring>
#include <vector>
#if defined(_WIN32)
#include <winsock2.h>
#include <ws2tcpip.h>
#endif
#if defined(__linux__)
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>
using SOCKET = int;
constexpr SOCKET INVALID_SOCKET = -1;
constexpr int SOCKET_ERROR = -1;
#endif

#if defined(_WIN32)
#pragma comment(lib, "ws2_32.lib")
#endif

namespace lzd_tools{
    class Socket {
    public:
        Socket(const std::string& ip, int port);
        ~Socket();

        bool connectToServer();
        bool SendDoubleArray(const std::vector<double>& data);

        /**
         * @brief 
         */
        bool ReceiveDoubleArray(std::vector<double>& data, size_t expectedSize);
        bool ReceiveIntArray(std::vector<int>& data, size_t expectedSize);

        bool SendJson(std::string json);
        bool ReceiveJson(std::string& json);

    private:
        SOCKET sock;
        struct sockaddr_in server_addr;

        bool initializeWinsock();
        void cleanupWinsock();
    };
}

