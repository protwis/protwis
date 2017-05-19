
from django.core.management.base import BaseCommand
import os

class Command(BaseCommand):

	def add_arguments(self, parser):
		parser.add_argument('stat_file', type=str, nargs='+')

	def handle(self, *args, **options):
		c = CountTemplate()
		c.count_rotamer_templates(options['stat_file'][0])

class CountTemplate():
	def __init__(self):
		pass

	def count_rotamer_templates(self, stat_file):
		templates = []
		with open(stat_file, 'r') as f:
			lines = f.readlines()
			for l in lines:
				line = l.split(',')
				if len(line)==5:
					temp = line[-1][13:17]
					if temp=='':
						continue
					if temp not in templates:
						print(line)
						templates.append(temp)
		print(templates)
		print(len(templates))
		return templates