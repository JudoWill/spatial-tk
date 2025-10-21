
venv:
	# Create a new conda environment in the venv directory
	conda create -p ./venv python=3.12
	conda run -p ./venv pip install -r requirements.txt

