#!/usr/bin/python

# this is my first go at doing stuff

import click, os, time, subprocess, shutil

@click.command()
@click.option('--count', default=1, help='Number of greetings.')
@click.option('--name', prompt='Your name',
              help='The person to greet.')
def hello(count, name):
    """Simple program that greets NAME for a total of COUNT times."""
    for x in range(count):
        click.echo('Hello '+ name + '!')

if __name__ == '__main__':
    hello()

print(cmd = 'ls')
