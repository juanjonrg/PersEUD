#!/usr/bin/env python3

import zmq, json, logging, time
from datetime import datetime
from random import randint

def dprint(*args, **kw):
    print('[' + datetime.now().isoformat(sep=' ', timespec='milliseconds') + ']', *args, **kw)

FRONTEND = "tcp://bullxual:3555"
REQUEST_TIMEOUT = 2000000
REQUEST_RETRIES = 100

dprint(f'Connecting to {FRONTEND}...')
context = zmq.Context()
client = context.socket(zmq.REQ)
client.setsockopt_string(zmq.IDENTITY, 'dummy{:04d}'.format(randint(0, 9999)))
client.connect(FRONTEND)

with open('params.json') as f:
    params = json.load(f)
client.send_string(json.dumps(params))

dprint('Parameters sent, waiting for a response...')
retries_left = REQUEST_RETRIES

while True:
    if (client.poll(REQUEST_TIMEOUT) & zmq.POLLIN) != 0:
        reply = client.recv()
        print("Server replied: {}".format(reply));
        retries_left = REQUEST_RETRIES
        break

    retries_left -= 1
    logging.warning("No response from server")
    # Socket is confused. Close and remove it.
    client.setsockopt(zmq.LINGER, 0)
    client.close()
    if retries_left == 0:
        logging.error("Server seems to be offline, abandoning")
        sys.exit()

    logging.info("Reconnecting to server…")
    # Create new connection
    client = context.socket(zmq.REQ)
    client.connect(SERVER_ENDPOINT)
    logging.info("Resending (%s)", request)
    client.send(request)

