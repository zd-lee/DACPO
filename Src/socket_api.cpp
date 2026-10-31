#include "socket_api.h"

#include <algorithm>
#include <cerrno>
#include <limits>
#include <nlohmann/json.hpp>

namespace {
int socket_last_error() {
#if defined(_WIN32)
    return WSAGetLastError();
#else
    return errno;
#endif
}

void close_socket(SOCKET sock) {
#if defined(_WIN32)
    closesocket(sock);
#else
    close(sock);
#endif
}

bool send_all(SOCKET sock, const char* data, size_t size) {
    size_t sent = 0;
    while (sent < size) {
        const int chunk = static_cast<int>((std::min)(size - sent,
            static_cast<size_t>((std::numeric_limits<int>::max)())));
#if defined(__linux__)
        const int count = send(sock, data + sent, chunk, MSG_NOSIGNAL);
#else
        const int count = send(sock, data + sent, chunk, 0);
#endif
        if (count <= 0) return false;
        sent += count;
    }
    return true;
}

bool receive_all(SOCKET sock, char* data, size_t size) {
    size_t received = 0;
    while (received < size) {
        const int chunk = static_cast<int>((std::min)(size - received,
            static_cast<size_t>((std::numeric_limits<int>::max)())));
        const int count = recv(sock, data + received, chunk, 0);
        if (count <= 0) return false;
        received += count;
    }
    return true;
}
}

lzd_tools::Socket::Socket(const std::string& ip, int port) : sock(INVALID_SOCKET) {
    if (!initializeWinsock()) {
        std::cerr << "Failed to initialize socket API." << std::endl;
        exit(EXIT_FAILURE);
    }

    sock = socket(AF_INET, SOCK_STREAM, 0);
    if (sock == INVALID_SOCKET) {
        std::cerr << "Socket creation error: " << socket_last_error() << std::endl;
        cleanupWinsock();
        exit(EXIT_FAILURE);
    }

    server_addr.sin_family = AF_INET;
    server_addr.sin_port = htons(port);
    if (inet_pton(AF_INET, ip.c_str(), &server_addr.sin_addr) <= 0) {
        std::cerr << "Invalid address/ Address not supported." << std::endl;
        cleanupWinsock();
        exit(EXIT_FAILURE);
    }
}

lzd_tools::Socket::~Socket() {
    if (sock != INVALID_SOCKET) {
        close_socket(sock);
    }
    cleanupWinsock();
}

bool lzd_tools::Socket::initializeWinsock() {
#if defined(_WIN32)
    WSADATA wsaData;
    return WSAStartup(MAKEWORD(2, 2), &wsaData) == 0;
#else
    return true;
#endif
}

void lzd_tools::Socket::cleanupWinsock() {
#if defined(_WIN32)
    WSACleanup();
#endif
}

bool lzd_tools::Socket::connectToServer() {
    if (connect(sock, (struct sockaddr*)&server_addr, sizeof(server_addr)) < 0) {
        std::cerr << "Connection Failed: " << socket_last_error() << std::endl;
        return false;
    }
    return true;
}

bool lzd_tools::Socket::SendDoubleArray(const std::vector<double>& data) {
    return send_all(sock, reinterpret_cast<const char*>(data.data()), data.size() * sizeof(double));
}

bool lzd_tools::Socket::ReceiveDoubleArray(std::vector<double>& data, size_t expectedSize) {
    data.resize(expectedSize);
    int bytesRead = recv(sock, reinterpret_cast<char*>(data.data()), expectedSize * sizeof(double), 0);
    if (bytesRead == SOCKET_ERROR) {
        std::cerr << "Receive failed: " << socket_last_error() << std::endl;
        return false;
    }
    return true;
}

bool lzd_tools::Socket::ReceiveIntArray(std::vector<int> &data, size_t expectedSize)
{
    data.resize(expectedSize);
    return receive_all(sock, reinterpret_cast<char*>(data.data()), expectedSize * sizeof(int));
}

bool lzd_tools::Socket::SendJson(std::string json)
{
    return send_all(sock, json.data(), json.size());
}

bool lzd_tools::Socket::ReceiveJson(std::string& json)
{
    json.clear();
    char buffer[1024];
    while (json.size() < 1024 * 1024) {
        const int count = recv(sock, buffer, sizeof(buffer), 0);
        if (count <= 0) return false;
        json.append(buffer, count);
        if (nlohmann::json::accept(json)) return true;
    }
    return false;
}
