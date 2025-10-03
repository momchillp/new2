
 // Initialize EmailJS
emailjs.init("NbSB7ErxhgVwuvkxk"); // Your public key

document.getElementById("chatbox-button").onclick = () => {
  document.getElementById("chatbox").style.display = "block";
};

document.getElementById("close-chat").onclick = () => {
  document.getElementById("chatbox").style.display = "none";
};

document.getElementById("send-message").onclick = () => {
  const name = document.getElementById("chat-name").value.trim();
  const email = document.getElementById("chat-email").value.trim();
  const message = document.getElementById("chat-message").value.trim();

  if (!name || !email || !message) {
    return alert("Please fill in all fields!");
  }

  emailjs.send("service_jxp4wim", "template_59tis3r", {
    name: name,
    email: email,
    message: message
  })
  .then(() => {
    alert("Message sent successfully!");
    document.getElementById("chat-name").value = "";
    document.getElementById("chat-email").value = "";
    document.getElementById("chat-message").value = "";
  })
  .catch((err) => {
    console.error("FAILED...", err);
    alert("Failed to send: " + JSON.stringify(err));
  });
};
