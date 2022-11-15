#!/usr/bin/env python3
# @Date       : 2022-10-12
# @Author     : Serafim Nenarokov (serafim.nenarokov@gmail.com)
# @Description: Simple wrapper around outbreak.info API for our purposes

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


class OutbreakApiError(Exception):
    def __init__(self, message):
        super().__init__(f"Outbreak API usage error: {message}")


class OutbreakUnauthorizedError(Exception):
    def __init__(self):
        super().__init__(f"Unauthorized")


class OutbreakApi:
    CURATED_LINEAGES_URL = "https://raw.githubusercontent.com/outbreak-info/outbreak.info/master/web/src/assets/genomics/curated_lineages.json"
    BASE_URL = "https://api.outbreak.info/genomics"
    AUTH_METHOD_NAME = "get-auth-token"


    def __init__(self, timeout=5, quiet=True):
        self.auth_key = None
        self.timeout = timeout
        self.quiet = quiet


    def auth_user(self):
        return self.__authorize()


    def prevalence(self, pangolin_lineage=None, location=None, mutations=None, cumulative=False):
        if pangolin_lineage == None and mutations == None:
            raise OutbreakApiError("Either pangolin_lineage or mutations should be specified.")

        cumulative = str(cumulative).lower()

        params = [['pangolin_lineage', pangolin_lineage],
                  ['location', location],
                  ['mutations', mutations],
                  ['cumulative', cumulative]]

        params = [p for p in params if p[1] != None]
        params = dict(params)

        return self.get_data('prevalence-by-location', params)


    def mutations_by_lineage(self, pangolin_lineage, frequency=0.75):
        params = {
            'pangolin_lineage': pangolin_lineage,
            'frequency': frequency
        }
        return self.get_data('lineage-mutations', params)


    def get_curated_lineages(self):
        url = OutbreakApi.CURATED_LINEAGES_URL
        return OutbreakApi.__raw_request(url, timeout=self.timeout, quiet=self.quiet)


    def get_data(self, method_name, params):
        if not self.__authorize():
            return None

        headers = {}
        headers['Authorization'] = f"Bearer {self.auth_key}"

        url = f"{self.BASE_URL}/{method_name}"

        result = OutbreakApi.__raw_request(url, params=params, headers=headers, timeout=self.timeout, quiet=self.quiet)

        if result and result['success'] == True:
            return result['results']


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
        response = OutbreakApi.__raw_request(auth_url, method='post', timeout=self.timeout, quiet=self.quiet)

        token = response['authn_token']
        url = response['authn_url']

        print("New API key is generated.")
        print("To activate it, please, open the following link, login to GISAID and close the browser window.")
        print(url)

        return token


    @staticmethod
    def __raw_request(url, method='get', params={}, headers={}, timeout=5, quiet=True):
        # standard params arg doesn't work with long values
        if len(params) != 0:
            params_str = '&'.join([f'{k}={v}' for k, v in params.items()])
            url = '?'.join([url, params_str])

        if not quiet:
            print(f"Requesting {method.upper()} `{url}`, timeout={timeout}, headers=`{headers}`")

        try:
            if method == 'get':
                response = requests.get(url, headers=headers, timeout=timeout)
            elif method == 'post':
                response = requests.post(url, headers=headers, timeout=timeout)
        except ConnectionError:
            raise OutbreakConnectionError()

        if not response.ok:
            raise OutbreakRequestError(f"{response.status_code} {response.reason}")

        return response.json()


