import socket

if __name__ == "__main__":
    ip = "localhost"
    port = 8888

    server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    server.bind((ip,port))
    server.listen(5)

    while True:
        client, address =  server.accept()
        print(f"Connection established - {address[0]}:{address[1]}")

        client.sendall("ACK!")
        #string = client.rec(1024)
        #string = string.decode("utf-8")
        #print(string)

        client.close()
