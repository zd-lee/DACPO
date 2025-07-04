#include "socket_api.h"

lzd_tools::Socket::Socket(const std::string& ip, int port) : sock(INVALID_SOCKET) {
    if (!initializeWinsock()) {
        std::cerr << "Failed to initialize Winsock." << std::endl;
        exit(EXIT_FAILURE);
    }

    sock = socket(AF_INET, SOCK_STREAM, 0);
    if (sock == INVALID_SOCKET) {
        std::cerr << "Socket creation error: " << WSAGetLastError() << std::endl;
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
    closesocket(sock);
    cleanupWinsock();
}

bool lzd_tools::Socket::initializeWinsock() {
    WSADATA wsaData;
    return WSAStartup(MAKEWORD(2, 2), &wsaData) == 0;
}

void lzd_tools::Socket::cleanupWinsock() {
    WSACleanup();
}

bool lzd_tools::Socket::connectToServer() {
    if (connect(sock, (struct sockaddr*)&server_addr, sizeof(server_addr)) < 0) {
        std::cerr << "Connection Failed: " << WSAGetLastError() << std::endl;
        return false;
    }
    return true;
}

bool lzd_tools::Socket::SendDoubleArray(const std::vector<double>& data) {
    int byteSent = 0;
    while(byteSent < data.size() * sizeof(double)) {
		int bytesSent = send(sock, reinterpret_cast<const char*>(data.data()) + byteSent, data.size() * sizeof(double) - byteSent, 0);
		if (bytesSent == SOCKET_ERROR) {
			std::cerr << "Send failed: " << WSAGetLastError() << std::endl;
			return false;
        }
        else {
			// std::cout << "Sent " << bytesSent << " bytes" << std::endl;
        }
		byteSent += bytesSent;
	}
    return true;
}

bool lzd_tools::Socket::ReceiveDoubleArray(std::vector<double>& data, size_t expectedSize) {
    data.resize(expectedSize);
    int bytesRead = recv(sock, reinterpret_cast<char*>(data.data()), expectedSize * sizeof(double), 0);
    if (bytesRead == SOCKET_ERROR) {
        std::cerr << "Receive failed: " << WSAGetLastError() << std::endl;
        return false;
    }
    return true;
}

bool lzd_tools::Socket::ReceiveIntArray(std::vector<int> &data, size_t expectedSize)
{
    data.resize(expectedSize);
    int bytesRead = recv(sock, reinterpret_cast<char*>(data.data()), expectedSize * sizeof(int), 0);
    if (bytesRead == SOCKET_ERROR) {
        std::cerr << "Receive failed: " << WSAGetLastError() << std::endl;
        return false;
    }
    return true;
}

bool lzd_tools::Socket::SendJson(std::string json)
{
    int bytesSent = send(sock, json.c_str(), json.size(), 0);
	if (bytesSent == SOCKET_ERROR) {
		std::cerr << "Send failed: " << WSAGetLastError() << std::endl;
		return false;
	}
	return true;
}

bool lzd_tools::Socket::ReceiveJson(std::string& json)
{
    json.clear();
    char buffer[1024];
	int bytesRead = recv(sock, buffer, 1024, 0);
	if (bytesRead == SOCKET_ERROR) {
		std::cerr << "Receive failed: " << WSAGetLastError() << std::endl;
		return false;
	}
	json = std::string(buffer, bytesRead);
	return true;
}