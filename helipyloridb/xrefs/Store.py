import json, requests
from enum import Enum


class Store:

    def __init__(self, serviceName):
        self.serviceName = serviceName

    def _check_provided(self, tocheck):

        for x in tocheck:
            if not x in self.provides:
                return False

        return True

    def _must_provide(self, tocheck):

        if not self._check_provided(tocheck):
            raise StoreException("Not all entities are provided by this store! \n Provided %s \n Requested %s" % (self.provides, tocheck))

        return True


class StoreException(Exception):

    def __init__(self, msg):
        self.msg = msg


class RequestMethod(Enum):

    GET = requests.get
    POST = requests.post

class RESTStore(Store):

    def __init__(self, serviceName, URL):

        super(RESTStore, self).__init__( serviceName )

        self.URL = URL

    def _requestGET(self, queryStr, params):

        resp = requests.get(url=self.URL + queryStr, params=params)
        data = resp.json()

        return data

    def _request(self, methodEnum, URLattach, args=None, json=None):

        response = methodEnum( self.URL + URLattach, args, json )

        return response

    def convert(self, fromEntity, toEntity, elements):
        pass

    def fetch(self, fromEntity, elements):
        pass
