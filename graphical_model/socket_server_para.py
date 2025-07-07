import socket
import numpy as np
import json
from pathlib import Path

# 设置服务器的IP和端口
HOST = '0.0.0.0'  # 监听所有IP地址
PORT = 11111    # 监听的端口号
REQUEST_BUFFER_SIZE = 1000 # 接收缓冲区大小，单位为字节
max_thread = 50 # 同时处理的最大线程数


import numpy as np
import gurobipy as gp
from tools import cal_loss

'''
data: 长度为n*n*2的数组
w为n*n的矩阵, w[i,j]表示x[i]和x[j]在指派相同的时候的权重,即A
inv_w为n*n的矩阵, inv_w[i,j]表示x[i]和x[j]在指派不同的时候的权重,即B
将w和inv_w平铺到一维数组后拼接,即data
'''
def data_preprocess(data):
    data = np.array(data)
    n = np.sqrt(data.shape[0] / 2).astype(np.int32)
    assert n * n * 2 == data.shape[0]
    data = data.reshape(2,n,n)
    w = data[0]
    inv_w = data[1]
    for i in range(n):
        assert w[i,i] == 0
        assert inv_w[i,i] == 0
        for j in range(n):
            # assert w[i,j] == w[j,i]
            # assert inv_w[i,j] == inv_w[j,i]
            if w[i,j] != 0:
                assert inv_w[i,j] != 0 or w[i,j] == 1
            elif inv_w[i,j] != 0:
                assert w[i,j] != 0 or inv_w[i,j] == 1
            if w[i,j] != 0 or inv_w[i,j] != 0:
                assert w[i,j] != inv_w[i,j]
    return w, inv_w



def MIQP(A,B):
    assert A.shape == B.shape
    assert A.shape[0] == A.shape[1]
    # Create a new model
    m = gp.Model("mip1")
    # Create variables
    n = len(A)
    x = m.addVars(n, vtype=gp.GRB.BINARY, name="x")
    # Set objective
    obj = gp.QuadExpr()
    obj -= cal_loss(x,A,B)
    m.setObjective(obj, gp.GRB.MAXIMIZE)

    # find the optimal solution
    m.optimize()
    res = np.zeros(n)
    
    # print('Obj: %g' % m.objVal)
    # for v in m.getVars():
    #     print('%s %g' % (v.varName, v.x))
    # print('Optimal solution found')
    
    print('Obj: %g' % m.objVal)
    for i in range(n):
        res[i] = x[i].x
    return res

def handle_client(conn, addr):
    with conn:
        print(f"Connected by {addr}")
        try:
            # 接收数据
            req = conn.recv(REQUEST_BUFFER_SIZE)
            req = json.loads(req.decode())
            print(req)
            data_buffer_size = req['data_size'] * 8
            repo = json.dumps({"status": "OK"})
            conn.sendall(repo.encode())
            data_recv = 0
            data = b''
            while data_recv < data_buffer_size:
                tdata = conn.recv(data_buffer_size - data_recv)
                data_recv += len(tdata)
                if not tdata:
                    break
                print(f"Received {len(tdata)} bytes")
                data += tdata
            if not data:
                return
            print(f"Received {len(data)} bytes in total")
            if len(data) != data_buffer_size:
                print(f"Data size mismatch. Expected {data_buffer_size} bytes, but received {len(data)} bytes.")
                assert False
    
            # 假设接收到的数据是一个numpy数组
            w, inv_w = data_preprocess(np.frombuffer(data, dtype=np.float64))
            res = MIQP(w, inv_w)
            # 返回结果
            conn.sendall(res.astype(np.int32).tobytes())
        except Exception as e:
            # 打印错误行号
            import traceback
            traceback.print_exc()
            print(f"Error: {e}")
            conn.sendall(json.dumps({"status": "ERROR"}).encode())
        finally:
            conn.close()
      
import threading
import time

def main():
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind((HOST, PORT))
        s.listen()
        print(f"Server listening on {HOST}:{PORT}")        
        while True:
            conn, addr = s.accept()
            while threading.active_count() > max_thread:
                time.sleep(1)
            t = threading.Thread(target=handle_client, args=(conn, addr))
            t.start()
            print(f"Active threads: {threading.active_count()}")
                 
if __name__ == "__main__":
    main()
