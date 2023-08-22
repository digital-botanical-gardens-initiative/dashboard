# Use an official Node.js runtime as the parent image
FROM node:14

# Set the working directory in the container
WORKDIR /usr/src/app

# Copy the current directory contents into the container at /usr/src/app
COPY package*.json ./

# Install the application's dependencies inside the container
RUN npm install

# If you want to install global packages, do it here
# For example: RUN npm install -g <package-name>

# Copy the rest of your app's source code from your host to your image filesystem.
COPY . .

# Make port 3000 available to the world outside this container
EXPOSE 3000

# Define the command to run your app using CMD which keeps the container running
CMD [ "node", "app.js" ]
