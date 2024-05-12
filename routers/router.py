#!/usr/bin/env python3

import zmq, json
from datetime import datetime

def dprint(*args, **kw):
    print('[' + datetime.now().isoformat(sep=' ', timespec='milliseconds') + ']', *args, **kw)

FRONTEND = "tcp://*:3559"
BACKEND = "tcp://*:3560"
context = zmq.Context(1)
frontend = context.socket(zmq.ROUTER)
backend = context.socket(zmq.ROUTER)
frontend.bind(FRONTEND) # For clients
backend.bind(BACKEND)  # For workers

poll = zmq.Poller()
poll.register(frontend, zmq.POLLIN)
poll.register(backend, zmq.POLLIN)

jobs = []
workers = []

dprint(f'Router started. Frontend {FRONTEND}, backend {BACKEND}.')

while True:
    dprint(f'Pending jobs: {len(jobs)}.')
    dprint(f'Available workers: {len(workers)} {workers}')

    socks = dict(poll.poll())
    # Worker management in the backend
    if socks.get(backend) == zmq.POLLIN:
        # Use worker address for LRU routing
        msg = backend.recv_multipart()
        if not msg:
            break
        address = msg[0]
        workers.append(address)

        # Everything after the second (delimiter) frame is reply
        reply = msg[2:]

        # Forward message to client if it's not a READY
        if reply[0] != b'READY':
            print("REPLY:", reply)
            frontend.send_multipart(reply)
            dprint('Worker {} replied: {}'.format(msg[0].decode(), msg[4].decode().rstrip()))
        else:
            dprint('Worker {} is ready.'.format(msg[0].decode()))

    # Client management in the frontend
    if socks.get(frontend) == zmq.POLLIN:
        msg = frontend.recv_multipart()
        print("FRONT END", msg)
        jobs.append(msg)

    # Job management
    if len(jobs) > 0 and len(workers) > 0:
        job = jobs.pop(0)
        worker = workers.pop(0)
        dprint(f'Assigning job from client {job[0]} to worker {worker.decode()}.')
        request = [worker, ''.encode()] + job
        backend.send_multipart(request)















