#!/usr/bin/env python3

"""
File:         download_Mathys2019_from_synapse.py
Created:      2022/09/21
Last Changed:
Author:       M.Vochteloo

Copyright (C) 2022 M.Vochteloo
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License can be found in the LICENSE file in the
root directory of this source tree. If not, see <https://www.gnu.org/licenses/>.
"""

# Standard imports.
import getpass

# Third party imports.
import synapseclient
import synapseutils

# Log in on Synapse.
user = input("Synapse username:")
password = getpass.getpass('Synapse password:')

syn = synapseclient.Synapse()
syn.login(user, password)

# Download.
_ = synapseutils.syncFromSynapse(syn,
                                 entity='syn18681734',
                                 path='processed/')
_ = synapseutils.syncFromSynapse(syn,
                                 entity='syn18638475',
                                 path='fastq/')
_ = synapseutils.syncFromSynapse(syn,
                                 entity='syn18642926',
                                 path='metadata/')