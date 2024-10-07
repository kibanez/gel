.PHONY: all setup install clean

all: setup install

setup:
		python3 -m venv .venv

install:
		.venv/bin/pip install -r requirements.txt
clean: 
		rm -rf .venv		