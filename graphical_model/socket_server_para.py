"""DACPO's dense socket protocol backed by OR-Tools CP-SAT."""

import argparse
import json
import math
import socketserver
import threading

import numpy as np

from cp_sat_optimizer import solve_flip, validate_costs


def matrix_dimension(size):
    if isinstance(size, bool) or not isinstance(size, int) or size < 0:
        raise ValueError("data_size must be a nonnegative integer")
    n = math.isqrt(size // 2)
    if 2 * n * n != size:
        raise ValueError("Expected data_size = 2 * n * n")
    return n


def data_preprocess(data):
    data = np.asarray(data, dtype=np.float64)
    n = matrix_dimension(data.size)
    return validate_costs(*data.reshape(2, n, n))


def recv_exact(conn, size):
    data = bytearray(size)
    view = memoryview(data)
    position = 0
    while position < size:
        received = conn.recv_into(view[position:])
        if received == 0:
            raise ConnectionError("Connection closed before receiving full payload")
        position += received
    return data


def recv_request(conn):
    data = bytearray()
    while len(data) < 4096:
        chunk = conn.recv(4096 - len(data))
        if not chunk:
            raise ConnectionError("Connection closed while reading JSON request")
        data.extend(chunk)
        try:
            request = json.loads(data)
        except (json.JSONDecodeError, UnicodeDecodeError):
            continue
        if not isinstance(request, dict):
            raise ValueError("Expected a JSON object")
        return request
    raise ValueError("Request header exceeds 4096 bytes")


class OptimizationHandler(socketserver.BaseRequestHandler):
    def handle(self):
        acknowledged = False
        self.request.settimeout(self.server.options.time_limit + 60)
        try:
            request = recv_request(self.request)
            if request.get("function_name") != "MIQP":
                raise ValueError("Only the MIQP function is supported")
            n = matrix_dimension(request.get("data_size"))
            if n > self.server.options.max_nodes:
                raise ValueError("Graph exceeds --max-nodes")
            self.request.sendall(b'{"status":"OK"}')
            acknowledged = True
            data = recv_exact(self.request, request["data_size"] * 8)
            same, opposite = data_preprocess(np.frombuffer(data, dtype="=f8"))
            with self.server.solve_slots:
                labels, stats = solve_flip(
                    same, opposite,
                    time_limit=self.server.options.time_limit,
                    workers=self.server.options.workers,
                    allow_feasible=self.server.options.allow_feasible,
                )
            print(json.dumps(stats), flush=True)
            self.request.sendall(labels.astype("=i4").tobytes())
        except Exception as exc:
            print(f"Optimization error from {self.client_address}: {exc}", flush=True)
            # After ACK the client expects binary labels: EOF signals failure.
            if not acknowledged:
                try:
                    self.request.sendall(json.dumps({"status": "ERROR", "error": str(exc)}).encode())
                except OSError:
                    pass


class OptimizationServer(socketserver.ThreadingTCPServer):
    allow_reuse_address = True
    daemon_threads = True

    def __init__(self, address, options):
        self.options = options
        self.solve_slots = threading.BoundedSemaphore(options.max_clients)
        super().__init__(address, OptimizationHandler)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=11111)
    parser.add_argument("--time-limit", type=float, default=300.0)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--max-clients", type=int, default=4)
    parser.add_argument("--max-nodes", type=int, default=4096)
    parser.add_argument("--allow-feasible", action="store_true")
    options = parser.parse_args()
    if not math.isfinite(options.time_limit) or options.time_limit <= 0 or not 1 <= options.workers <= 32:
        parser.error("time-limit must be positive and workers must be in [1, 32]")
    if options.max_clients < 1 or options.max_nodes < 1:
        parser.error("max-clients and max-nodes must be positive")
    with OptimizationServer((options.host, options.port), options) as server:
        print(f"DACPO CP-SAT listening on {options.host}:{server.server_address[1]}", flush=True)
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            pass


if __name__ == "__main__":
    main()
