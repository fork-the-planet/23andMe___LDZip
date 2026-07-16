# GitHub Actions Workflows

This directory contains the GitHub Actions workflows used by the LDZip repository.

## Workflows

### `ci.yaml`

Runs the LDZip automated tests on:

- Pushes to `main`
- Pull requests targeting `main`

All jobs run inside the shared Docker image:

```text
ghcr.io/23andme/ldzip/ldzip-ci:latest
```

### `build-image.yaml`

Builds and publishes the CI Docker image to GitHub Container Registry (GHCR):

```text
ghcr.io/23andme/ldzip/ldzip-ci:latest
```

This workflow is currently disabled (via `if: false` in the `build-and-push` job) because the GitHub organization does not allow the third-party Docker GitHub Actions it requires. Until the workflow is enabled, the CI Docker image must be built and published manually.

## Manually Updating the CI Docker Image

Run the following commands from the root of the LDZip repository (where `Dockerfile.ci` is located):

```bash
docker login ghcr.io -u USERNAME

docker buildx build \
  --platform linux/amd64 \
  -f Dockerfile.ci \
  -t ghcr.io/23andme/ldzip/ldzip-ci:latest \
  .

docker images | grep ldzip-ci

docker push ghcr.io/23andme/ldzip/ldzip-ci:latest
```

### Notes

- Replace `USERNAME` with your GitHub username.
- Authenticate using a GitHub Personal Access Token when prompted for your password.
- The token must have permission to publish packages to the `23andMe` organization.
- After the image is pushed, subsequent runs of `ci.yaml` will automatically use the updated image.