ifeq ($(OS), Windows_NT)
	pwd_path := %cd%
else
	pwd_path := ${PWD}
endif

build: 
	docker build -t cec2020 .

run:
	@echo $(pwd_path)
	docker run -it --rm -v $(pwd_path):/cec2020 cec2020 $$cec	

clean: 
	rm -r $(pwd_path)/data/reproduction
