from unittest import TestCase
from helpers.helpers import get_git_revision_hash

__author__ = 'hayssam'


class TestGet_git_revision_hash(TestCase):
	def test_get_git_revision_hash(self):
		print get_git_revision_hash()
