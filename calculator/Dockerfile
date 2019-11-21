FROM mcr.microsoft.com/dotnet/core/sdk:3.0

# Install prerequisites for System.Drawing.Common
RUN apt-get update && \
    apt-get install -y libc6-dev libgdiplus

# Environment variables
ENV test 1

# Install the application into /app
# COPY /home/mshirley/git/github/parallel_lighting/ /app

VOLUME ["/data", "/app"]

ENTRYPOINT dotnet run --project calculator
