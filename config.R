#### Mail Config
defaultUsers <- c(
	'<name@domain.com>'
	) # Add as necessary

email_list_for_structure <- list()
email_list_for_structure[['protein']] <- c(defaultUsers, '<anothername@anotherdomain.com>') # And so on (requires explicit knowledge of the protein list)
