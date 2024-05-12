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
submitted_jobs = 0

jobWorker = []
index_job = 0

dprint(f'Router started. Frontend {FRONTEND}, backend {BACKEND}.')

while True:
    print("a")
    dprint(f'Pending jobs: {len(jobs)}.')
    dprint(f'Available workers: {len(workers)} {workers}')

    
    socks = dict(poll.poll(5000))
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
                
            submitted_jobs -= 1
            print("MSG", msg)
            for i in range(0, len(jobWorker)):
                print("Entra", jobWorker[i][0], address)
                if jobWorker[i]==(address.decode("utf-8")):
                    print("Es igual")
                    jobWorker[i] = msg[2:]
        
            if submitted_jobs == 0:
                merged_reply = jobWorker[0]
                print("MERGED REPLY", merged_reply)
                print(len(jobWorker))

                for i in range(1,len(jobWorker)):

                    # Parse JSON strings to Python objects
                    obj_list1 = json.loads(merged_reply[2]) + json.loads(jobWorker[i][2])

                    # Convert the result to JSON string
                    merged_reply[2] = json.dumps(obj_list1).encode()
                
                frontend.send_multipart(merged_reply)
                dprint('Worker {} replied: {}'.format(msg[0].decode(), msg[4].decode().rstrip()))
                merged_reply = []

        else:
            dprint('Worker {} is ready.'.format(msg[0].decode()))

    # Client management in the frontend
    if socks.get(frontend) == zmq.POLLIN:
        msg = frontend.recv_multipart()
        
        # Decode bytes to string and load JSON
        json_obj = json.loads(msg[2].decode('utf-8'))

        # Specify the number of splits
        num_splits = 2  # Split into 2 parts, change this as needed

        # Calculate the step size for splitting
        num_plans = len(json_obj["Plans"])
        step = num_splits# num_plans // num_splits
        submitted_jobs = num_plans // num_splits#num_splits
        # Splitting the "Plans" list based on the step size
        plans_list = json_obj["Plans"]
        split_plans = [plans_list[i:i + step] for i in range(0, num_plans, step)]

        # Create dictionaries for each split
        split_json_objs = []
        for i, plans in enumerate(split_plans):
            split_obj = json_obj.copy()
            split_obj["Plans"] = plans
            split_json_objs.append(split_obj)

        # Convert dictionaries back to JSON strings
        json_strs = [json.dumps(obj,indent=None) for obj in split_json_objs]

        # Print or use the JSON strings as needed
        for i, json_str in enumerate(json_strs):
            msg[2]= json_str.encode()
            jobs.append(msg.copy())
            print("JOB APPEND", msg)
            print("JOBS AFTER APPEND", jobs)
            jobWorker.append("-")
            

    # Job management
    if len(jobs) > 0 and len(workers) > 0:
        job = jobs.pop(0)
        print("JOB:", job)
        print("JOBS AFTER POP", jobs)
        worker = workers.pop(0)
        jobWorker[index_job] = worker.decode("utf-8")
        index_job += 1
        dprint(f'Assigning job from client {job[0]} to worker {worker.decode()}.')
        request = [worker, ''.encode()] + job
        backend.send_multipart(request)
        















