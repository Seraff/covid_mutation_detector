#!/usr/bin/env python3
# @Date    : 2022-10-12
# @Author  : Serafim Nenarokov (serafim.nenarokov@gmail.com)

import os
import pathlib
import requests

ROOT_PATH = str(pathlib.Path(os.path.dirname(os.path.realpath(__file__))).parent)

AUTH_FILENAME = ".outbreak_api_key"
AUTH_FILE_PATH = os.path.join(ROOT_PATH, AUTH_FILENAME)


class OutbreakConnectionError(Exception):
    def __init__(self):
        super().__init__("Cannot connect to Outbreak API")


class OutbreakRequestError(Exception):
    def __init__(self, message):
        super().__init__(f"Outbreak API returns a bad status: {message}")


class OutbreakApi:
    BASE_URL = "https://api.outbreak.info/genomics"
    AUTH_METHOD_NAME = "get-auth-token"


    def __init__(self):
        self.auth_key = None


    def get_data(self, method_name, params):
        if not self.__authorize():
            return None

        headers = {}
        headers['Authorization'] = f"Bearer {self.auth_key}"

        url = f"{self.BASE_URL}/{method_name}"

        return OutbreakApi.__raw_request(url, params=params, headers=headers)


    def __authorize(self):
        """Make sure, that we are authorized"""

        if self.auth_key == None:
            if os.path.exists(AUTH_FILE_PATH):
                with open(AUTH_FILE_PATH) as f:
                    self.auth_key = f.read().strip()
                return True
            else:
                self.auth_key = self.__create_auth()

                with open(AUTH_FILE_PATH, 'w') as f:
                    f.write(self.auth_key)
                return False
        else:
            return True


    def __create_auth(self):
        """Request new bearer token"""

        auth_url = f"{self.BASE_URL}/{self.AUTH_METHOD_NAME}"
        response = OutbreakApi.__raw_request(auth_url, method='post', timeout=5)

        token = response['authn_token']
        url = response['authn_url']

        print("New API key is generated.")
        print("To activate it, please, open the following link, login to GISAID and close the browser window.")
        print(url)

        return token


    @staticmethod
    def __raw_request(url, method='get', params={}, headers={}, timeout=5):
        try:
            if method == 'get':
                response = requests.get(url, params=params, headers=headers, timeout=timeout)
            elif method == 'post':
                response = requests.post(url, params=params, headers=headers, timeout=timeout)
        except ConnectionError:
            raise OutbreakConnectionError()

        if not response.ok:
            raise OutbreakRequestError(f"{response.status_code} {response.reason}")

        return response.json()


