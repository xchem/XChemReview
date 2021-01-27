# Example config file.

#### DBConfig
db_name <- 'test'
host_db <- '0.0.0.0'
db_port <- '0000'
db_user <- 'username'
db_password <- 'password'

#### Slack Config
apiuser <- 'SlackUserAPIKey'
api <- 'SlackChannelAPIKey'

#### Mail Config
defaultUsers <- c(
	'<name@domain.com>'
	) # Add as necessary

email_list_for_structure <- list()
email_list_for_structure[['protein']] <- c(defaultUsers, '<anothername@anotherdomain.com>') # And so on (requires explicit knowledge of the protein list)